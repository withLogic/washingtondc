/*******************************************************************************
 *
 *
 *    WashingtonDC Dreamcast Emulator
 *    Copyright (C) 2017 snickerbockers
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 ******************************************************************************/

#ifndef SH4_INST_H_
#define SH4_INST_H_

#include <assert.h>
#include <stdint.h>
#include <stdbool.h>

#include "error.h"
#include "types.h"

struct Sh4;
typedef struct Sh4 Sh4;

union Sh4OpArgs {
    inst_t inst;

    // operators that take a single general-purpose register and/or
    // an 8-bit immediate value
    struct {
        inst_t imm8 : 8;
        inst_t gen_reg : 4;
    inst_t : 4;
    };

    // signed 8-bit immediate value
    struct {
        int8_t simm8 : 8;
    inst_t : 8;
    };

    // operators that take a base register and a 4-bit displacement
    struct {
        inst_t imm4 : 4;

        /* the nomenclature here is kind of confusing:
         * base_reg_src can be a source *or* a dest for opcodes that only
         * take one register.  For opcodes that need two registers,
         * base_reg_src is the source and base_reg_dst is the destination.
         */
        inst_t base_reg_src : 4;
        inst_t base_reg_dst : 4;
    };

    // operators that take a 12-bit immediate value
    struct {
        inst_t imm12 : 12;
    inst_t : 4;
    };

    // signed 12-bit immediate value
    struct {
        int16_t simm12 : 12;
    inst_t : 4;
    };

    // operators that take two general-purpose registers
    struct {
    inst_t : 4;
        inst_t src_reg : 4;
        inst_t dst_reg : 4;
    inst_t : 4;
    };

    // opcodes that take a single floating-point register
    struct {
    inst_t : 8;
        inst_t fr_reg : 4;
    inst_t : 4;
    };

    // opcodes that take two floating-point registers
    struct {
    inst_t : 4;
        inst_t fr_src : 4;
        inst_t fr_dst : 4;
    inst_t : 4;
    };

    // opcodes that take a single double-precision floating-point register
    struct {
    inst_t : 9;
        inst_t dr_reg : 3;
    inst_t : 4;
    };

    // opcodes that take two double-precision floating-point registers
    struct {
    inst_t : 5;
        inst_t dr_src : 3;
    inst_t : 1;
        inst_t dr_dst : 3;
    inst_t : 4;
    };

    // opcodes that take a single double-precision floating-point register
    struct {
    inst_t : 9;
        inst_t xd_reg : 3;
    inst_t : 4;
    };

    // opcodes that take two double-precision floating-point registers
    struct {
    inst_t : 5;
        inst_t xd_src : 3;
    inst_t : 1;
        inst_t xd_dst : 3;
    inst_t : 4;
    };

    // opcodes that take two floating-point vector registers
    struct {
    inst_t : 8;
        inst_t fv_src : 2;
        inst_t fv_dst : 2;
    inst_t : 4;
    };

    // opcodes that take a single floating-point vector register
    struct {
    inst_t : 10;
        inst_t fv_reg : 2;
    inst_t : 4;
    };

    // operators that take in a banked register
    struct {
    inst_t : 4;
        inst_t bank_reg : 3;
    inst_t : 9;
    };
};

typedef union Sh4OpArgs Sh4OpArgs;

/*
 * the lut is a static (global) table that will be shared by all sh4
 * instances even if there's more than one of them, but sh4_init_lut
 * will always initialize it in the exact same way so it's safe to call this
 * function more than once.
 */
void sh4_init_inst_lut();

static_assert(sizeof(Sh4OpArgs) == 2,
              "sizeof(Sh4OpArgs) must match the size of an sh4 instruction!");

typedef void (*opcode_func_t)(Sh4*, Sh4OpArgs oa);

/*
 * The Hitach SH-4 is a dual-issue CPU, meaning that there are two separate
 * pipelines capable of executing instructions simultaneously.  From the
 * software's perspective, instruction execution is sequential, so normal
 * pipeline limitations such as stalls can still apply.
 *
 * Assuming that there are no stalls, the rule is that there are 6 distinct
 * groups of instructions (see the group member of InstOpcode) and that what
 * group an instruction is in determines which other groups it can execute in
 * parallel with.  The MT group can execute in parallel with any instruction
 * group except for CO (even itself), CO cannot execute in parallel with any
 * group (not even itself) and all every group is capable of executing in
 * parallel with any group except for itself and CO.
 *
 * OBSERVATION:
 * every instruction that takes more than 1 cycle to execute is part of group
 * CO.  CO instructions never execute in parallel.  This makes the
 * cycle-counting significantly simpler because I know that I will never need to
 * model a situation where one of the pipelines is executing an instruction that takes
 * long that what the other pipeline is executing.
 */
typedef enum sh4_inst_group {
    SH4_GROUP_MT,
    SH4_GROUP_EX,
    SH4_GROUP_BR,
    SH4_GROUP_LS,
    SH4_GROUP_FE,
    SH4_GROUP_CO,

    /*
     * used by the sh4_single_step code to indicate that the previous
     * instruction was an "even" instruction, meaning that this instruction
     * will not be free under any circumstance (althoutgh the next one might).
     *
     * Obviously this is not a real instruction group.
     */
    SH4_GROUP_NONE
} sh4_inst_group_t;

struct InstOpcode {
    // format string compiled to make mask and val
    char const *fmt;

    // opcode handler function
    opcode_func_t func;

    // if this is true, this inst cant be called from a delay slot
    bool is_branch;

    /*
     * execution group.  If I was emulating the dual-issue nature of the
     * pipeline, this would determine which instruction could execute
     * simoltaneously
     */
    sh4_inst_group_t group;

    /*
     * Number of cycles after each instruction before the next instruction can
     * be issued within the same pipeline.  The other constraining factor that
     * can delay is the latency (how long it takes an instruction's output to
     * become available), but I don't store that because some opcodes don't have
     * uniform latency, and some opcodes have multiple latencies for different
     * outputs
     */
    unsigned issue;

    // instructions are matched to this opcode
    // by anding with mask and checking for equality with val
    inst_t mask;
    inst_t val;
};

typedef struct InstOpcode InstOpcode;

/*
 * maps 16-bit instructions to InstOpcodes for O(1) decoding
 * this array looks big but it's really only half a megabyte
 */
extern InstOpcode const *sh4_inst_lut[1 << 16];

void sh4_compile_instructions(Sh4 *sh4);
void sh4_compile_instruction(Sh4 *sh4, struct InstOpcode *op);

ERROR_STRING_ATTR(opcode_format);
ERROR_STRING_ATTR(opcode_name);
ERROR_INT_ATTR(instruction);
ERROR_INT_ATTR(instruction_mask);
ERROR_INT_ATTR(instruction_expect);
ERROR_U32_ATTR(fpscr);
ERROR_U32_ATTR(fpscr_expect);
ERROR_U32_ATTR(fpscr_mask);
ERROR_INT_ATTR(inst_bin);

#define SH4_INST_RAISE_ERROR(sh4, error_tp)     \
    do {                                        \
        RAISE_ERROR(error_tp);                  \
    } while (0)

#endif
