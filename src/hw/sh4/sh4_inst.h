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
extern opcode_func_t sh4_inst_lut[1 << 16];

void sh4_compile_instructions(Sh4 *sh4);
void sh4_compile_instruction(Sh4 *sh4, struct InstOpcode *op);

// RTS
// 0000000000001011
void sh4_inst_rts(Sh4 *sh4, Sh4OpArgs inst);

// CLRMAC
// 0000000000101000
void sh4_inst_clrmac(Sh4 *sh4, Sh4OpArgs inst);

// CLRS
// 0000000001001000
void sh4_inst_clrs(Sh4 *sh4, Sh4OpArgs inst);

// CLRT
// 0000000000001000
void sh4_inst_clrt(Sh4 *sh4, Sh4OpArgs inst);

// LDTLB
// 0000000000111000
void sh4_inst_ldtlb(Sh4 *sh4, Sh4OpArgs inst);

// NOP
// 0000000000001001
void sh4_inst_nop(Sh4 *sh4, Sh4OpArgs inst);

// RTE
// 0000000000101011
void sh4_inst_rte(Sh4 *sh4, Sh4OpArgs inst);

// SETS
// 0000000001011000
void sh4_inst_sets(Sh4 *sh4, Sh4OpArgs inst);

// SETT
// 0000000000011000
void sh4_inst_sett(Sh4 *sh4, Sh4OpArgs inst);

// SLEEP
// 0000000000011011
void sh4_inst_sleep(Sh4 *sh4, Sh4OpArgs inst);

// FRCHG
// 1111101111111101
void sh4_inst_frchg(Sh4 *sh4, Sh4OpArgs inst);

// FSCHG
// 1111001111111101
void sh4_inst_fschg(Sh4 *sh4, Sh4OpArgs inst);

// MOVT Rn
// 0000nnnn00101001
void sh4_inst_unary_movt_gen(Sh4 *sh4, Sh4OpArgs inst);

// CMP/PZ Rn
// 0100nnnn00010001
void sh4_inst_unary_cmppz_gen(Sh4 *sh4, Sh4OpArgs inst);

// CMP/PL Rn
// 0100nnnn00010101
void sh4_inst_unary_cmppl_gen(Sh4 *sh4, Sh4OpArgs inst);

// DT Rn
// 0100nnnn00010000
void sh4_inst_unary_dt_gen(Sh4 *sh4, Sh4OpArgs inst);

// ROTL Rn
// 0100nnnn00000100
void sh4_inst_unary_rotl_gen(Sh4 *sh4, Sh4OpArgs inst);

// ROTR Rn
// 0100nnnn00000101
void sh4_inst_unary_rotr_gen(Sh4 *sh4, Sh4OpArgs inst);

// ROTCL Rn
// 0100nnnn00100100
void sh4_inst_unary_rotcl_gen(Sh4 *sh4, Sh4OpArgs inst);

// ROTCR Rn
// 0100nnnn00100101
void sh4_inst_unary_rotcr_gen(Sh4 *sh4, Sh4OpArgs inst);

// SHAL Rn
// 0100nnnn00200000
void sh4_inst_unary_shal_gen(Sh4 *sh4, Sh4OpArgs inst);

// SHAR Rn
// 0100nnnn00100001
void sh4_inst_unary_shar_gen(Sh4 *sh4, Sh4OpArgs inst);

// SHLL Rn
// 0100nnnn00000000
void sh4_inst_unary_shll_gen(Sh4 *sh4, Sh4OpArgs inst);

// SHLR Rn
// 0100nnnn00000001
void sh4_inst_unary_shlr_gen(Sh4 *sh4, Sh4OpArgs inst);

// SHLL2 Rn
// 0100nnnn00001000
void sh4_inst_unary_shll2_gen(Sh4 *sh4, Sh4OpArgs inst);

// SHLR2 Rn
// 0100nnnn00001001
void sh4_inst_unary_shlr2_gen(Sh4 *sh4, Sh4OpArgs inst);

// SHLL8 Rn
// 0100nnnn00011000
void sh4_inst_unary_shll8_gen(Sh4 *sh4, Sh4OpArgs inst);

// SHLR8 Rn
// 0100nnnn00011001
void sh4_inst_unary_shlr8_gen(Sh4 *sh4, Sh4OpArgs inst);

// SHLL16 Rn
// 0100nnnn00101000
void sh4_inst_unary_shll16_gen(Sh4 *sh4, Sh4OpArgs inst);

// SHLR16 Rn
// 0100nnnn00101001
void sh4_inst_unary_shlr16_gen(Sh4 *sh4, Sh4OpArgs inst);

// BRAF Rn
// 0000nnnn00100011
void sh4_inst_unary_braf_gen(Sh4 *sh4, Sh4OpArgs inst);

// BSRF Rn
// 0000nnnn00000011
void sh4_inst_unary_bsrf_gen(Sh4 *sh4, Sh4OpArgs inst);

// CMP/EQ #imm, R0
// 10001000iiiiiiii
void sh4_inst_binary_cmpeq_imm_r0(Sh4 *sh4, Sh4OpArgs inst);

// AND.B #imm, @(R0, GBR)
// 11001101iiiiiiii
void sh4_inst_binary_andb_imm_r0_gbr(Sh4 *sh4, Sh4OpArgs inst);

// AND #imm, R0
// 11001001iiiiiiii
void sh4_inst_binary_and_imm_r0(Sh4 *sh4, Sh4OpArgs inst);

// OR.B #imm, @(R0, GBR)
// 11001111iiiiiiii
void sh4_inst_binary_orb_imm_r0_gbr(Sh4 *sh4, Sh4OpArgs inst);

// OR #imm, R0
// 11001011iiiiiiii
void sh4_inst_binary_or_imm_r0(Sh4 *sh4, Sh4OpArgs inst);

// TST #imm, R0
// 11001000iiiiiiii
void sh4_inst_binary_tst_imm_r0(Sh4 *sh4, Sh4OpArgs inst);

// TST.B #imm, @(R0, GBR)
// 11001100iiiiiiii
void sh4_inst_binary_tstb_imm_r0_gbr(Sh4 *sh4, Sh4OpArgs inst);

// XOR #imm, R0
// 11001010iiiiiiii
void sh4_inst_binary_xor_imm_r0(Sh4 *sh4, Sh4OpArgs inst);

// XOR.B #imm, @(R0, GBR)
// 11001110iiiiiiii
void sh4_inst_binary_xorb_imm_r0_gbr(Sh4 *sh4, Sh4OpArgs inst);

// BF label
// 10001011dddddddd
void sh4_inst_unary_bf_disp(Sh4 *sh4, Sh4OpArgs inst);

// BF/S label
// 10001111dddddddd
void sh4_inst_unary_bfs_disp(Sh4 *sh4, Sh4OpArgs inst);

// BT label
// 10001001dddddddd
void sh4_inst_unary_bt_disp(Sh4 *sh4, Sh4OpArgs inst);

// BT/S label
// 10001101dddddddd
void sh4_inst_unary_bts_disp(Sh4 *sh4, Sh4OpArgs inst);

// BRA label
// 1010dddddddddddd
void sh4_inst_unary_bra_disp(Sh4 *sh4, Sh4OpArgs inst);

// BSR label
// 1011dddddddddddd
void sh4_inst_unary_bsr_disp(Sh4 *sh4, Sh4OpArgs inst);

// TRAPA #immed
// 11000011iiiiiiii
void sh4_inst_unary_trapa_disp(Sh4 *sh4, Sh4OpArgs inst);

// TAS.B @Rn
// 0100nnnn00011011
void sh4_inst_unary_tasb_gen(Sh4 *sh4, Sh4OpArgs inst);

// OCBI @Rn
// 0000nnnn10010011
void sh4_inst_unary_ocbi_indgen(Sh4 *sh4, Sh4OpArgs inst);

// OCBP @Rn
// 0000nnnn10100011
void sh4_inst_unary_ocbp_indgen(Sh4 *sh4, Sh4OpArgs inst);

// OCBWB @Rn
// 0000nnnn10110011
void sh4_inst_unary_ocbwb_indgen(Sh4 *sh4, Sh4OpArgs inst);

// PREF @Rn
// 0000nnnn10000011
void sh4_inst_unary_pref_indgen(Sh4 *sh4, Sh4OpArgs inst);

// JMP @Rn
// 0100nnnn00101011
void sh4_inst_unary_jmp_indgen(Sh4 *sh4, Sh4OpArgs inst);

// JSR @Rn
//0100nnnn00001011
void sh4_inst_unary_jsr_indgen(Sh4 *sh4, Sh4OpArgs inst);

// LDC Rm, SR
// 0100mmmm00001110
void sh4_inst_binary_ldc_gen_sr(Sh4 *sh4, Sh4OpArgs inst);

// LDC Rm, GBR
// 0100mmmm00011110
void sh4_inst_binary_ldc_gen_gbr(Sh4 *sh4, Sh4OpArgs inst);

// LDC Rm, VBR
// 0100mmmm00101110
void sh4_inst_binary_ldc_gen_vbr(Sh4 *sh4, Sh4OpArgs inst);

// LDC Rm, SSR
// 0100mmmm00111110
void sh4_inst_binary_ldc_gen_ssr(Sh4 *sh4, Sh4OpArgs inst);

// LDC Rm, SPC
// 0100mmmm01001110
void sh4_inst_binary_ldc_gen_spc(Sh4 *sh4, Sh4OpArgs inst);

// LDC Rm, DBR
// 0100mmmm11111010
void sh4_inst_binary_ldc_gen_dbr(Sh4 *sh4, Sh4OpArgs inst);

// STC SR, Rn
// 0000nnnn00000010
void sh4_inst_binary_stc_sr_gen(Sh4 *sh4, Sh4OpArgs inst);

// STC GBR, Rn
// 0000nnnn00010010
void sh4_inst_binary_stc_gbr_gen(Sh4 *sh4, Sh4OpArgs inst);

// STC VBR, Rn
// 0000nnnn00100010
void sh4_inst_binary_stc_vbr_gen(Sh4 *sh4, Sh4OpArgs inst);

// STC SSR, Rn
// 0000nnnn01000010
void sh4_inst_binary_stc_ssr_gen(Sh4 *sh4, Sh4OpArgs inst);

// STC SPC, Rn
// 0000nnnn01000010
void sh4_inst_binary_stc_spc_gen(Sh4 *sh4, Sh4OpArgs inst);

// STC SGR, Rn
// 0000nnnn00111010
void sh4_inst_binary_stc_sgr_gen(Sh4 *sh4, Sh4OpArgs inst);

// STC DBR, Rn
// 0000nnnn11111010
void sh4_inst_binary_stc_dbr_gen(Sh4 *sh4, Sh4OpArgs inst);

// LDC.L @Rm+, SR
// 0100mmmm00000111
void sh4_inst_binary_ldcl_indgeninc_sr(Sh4 *sh4, Sh4OpArgs inst);

// LDC.L @Rm+, GBR
// 0100mmmm00010111
void sh4_inst_binary_ldcl_indgeninc_gbr(Sh4 *sh4, Sh4OpArgs inst);

// LDC.L @Rm+, VBR
// 0100mmmm00100111
void sh4_inst_binary_ldcl_indgeninc_vbr(Sh4 *sh4, Sh4OpArgs inst);

// LDC.L @Rm+, SSR
// 0100mmmm00110111
void sh4_inst_binary_ldcl_indgenic_ssr(Sh4 *sh4, Sh4OpArgs inst);

// LDC.L @Rm+, SPC
// 0100mmmm01000111
void sh4_inst_binary_ldcl_indgeninc_spc(Sh4 *sh4, Sh4OpArgs inst);

// LDC.L @Rm+, DBR
// 0100mmmm11110110
void sh4_inst_binary_ldcl_indgeninc_dbr(Sh4 *sh4, Sh4OpArgs inst);

// STC.L SR, @-Rn
// 0100nnnn00000011
void sh4_inst_binary_stcl_sr_inddecgen(Sh4 *sh4, Sh4OpArgs inst);

// STC.L GBR, @-Rn
// 0100nnnn00010011
void sh4_inst_binary_stcl_gbr_inddecgen(Sh4 *sh4, Sh4OpArgs inst);

// STC.L VBR, @-Rn
// 0100nnnn00100011
void sh4_inst_binary_stcl_vbr_inddecgen(Sh4 *sh4, Sh4OpArgs inst);

// STC.L SSR, @-Rn
// 0100nnnn00110011
void sh4_inst_binary_stcl_ssr_inddecgen(Sh4 *sh4, Sh4OpArgs inst);

// STC.L SPC, @-Rn
// 0100nnnn01000011
void sh4_inst_binary_stcl_spc_inddecgen(Sh4 *sh4, Sh4OpArgs inst);

// STC.L SGR, @-Rn
// 0100nnnn00110010
void sh4_inst_binary_stcl_sgr_inddecgen(Sh4 *sh4, Sh4OpArgs inst);

// STC.L DBR, @-Rn
// 0100nnnn11110010
void sh4_inst_binary_stcl_dbr_inddecgen(Sh4 *sh4, Sh4OpArgs inst);

// MOV #imm, Rn
// 1110nnnniiiiiiii
void sh4_inst_binary_mov_imm_gen(Sh4 *sh4, Sh4OpArgs inst);

// ADD #imm, Rn
// 0111nnnniiiiiiii
void sh4_inst_binary_add_imm_gen(Sh4 *sh4, Sh4OpArgs inst);

// MOV.W @(disp, PC), Rn
// 1001nnnndddddddd
void sh4_inst_binary_movw_binind_disp_pc_gen(Sh4 *sh4, Sh4OpArgs inst);

// MOV.L @(disp, PC), Rn
// 1101nnnndddddddd
void sh4_inst_binary_movl_binind_disp_pc_gen(Sh4 *sh4, Sh4OpArgs inst);

// MOV Rm, Rn
// 0110nnnnmmmm0011
void sh4_inst_binary_movw_gen_gen(Sh4 *sh4, Sh4OpArgs inst);

// SWAP.B Rm, Rn
// 0110nnnnmmmm1000
void sh4_inst_binary_swapb_gen_gen(Sh4 *sh4, Sh4OpArgs inst);

// SWAP.W Rm, Rn
// 0110nnnnmmmm1001
void sh4_inst_binary_swapw_gen_gen(Sh4 *sh4, Sh4OpArgs inst);

// XTRCT Rm, Rn
// 0110nnnnmmmm1101
void sh4_inst_binary_xtrct_gen_gen(Sh4 *sh4, Sh4OpArgs inst);

// ADD Rm, Rn
// 0011nnnnmmmm1100
void sh4_inst_binary_add_gen_gen(Sh4 *sh4, Sh4OpArgs inst);

// ADDC Rm, Rn
// 0011nnnnmmmm1110
void sh4_inst_binary_addc_gen_gen(Sh4 *sh4, Sh4OpArgs inst);

// ADDV Rm, Rn
// 0011nnnnmmmm1111
void sh4_inst_binary_addv_gen_gen(Sh4 *sh4, Sh4OpArgs inst);

// CMP/EQ Rm, Rn
// 0011nnnnmmmm0000
void sh4_inst_binary_cmpeq_gen_gen(Sh4 *sh4, Sh4OpArgs inst);

// CMP/HS Rm, Rn
// 0011nnnnmmmm0010
void sh4_inst_binary_cmphs_gen_gen(Sh4 *sh4, Sh4OpArgs inst);

// CMP/GE Rm, Rn
// 0011nnnnmmmm0011
void sh4_inst_binary_cmpge_gen_gen(Sh4 *sh4, Sh4OpArgs inst);

// CMP/HI Rm, Rn
// 0011nnnnmmmm0110
void sh4_inst_binary_cmphi_gen_gen(Sh4 *sh4, Sh4OpArgs inst);

// CMP/GT Rm, Rn
// 0011nnnnmmmm0111
void sh4_inst_binary_cmpgt_gen_gen(Sh4 *sh4, Sh4OpArgs inst);

// CMP/STR Rm, Rn
// 0010nnnnmmmm1100
void sh4_inst_binary_cmpstr_gen_gen(Sh4 *sh4, Sh4OpArgs inst);

// DIV1 Rm, Rn
// 0011nnnnmmmm0100
void sh4_inst_binary_div1_gen_gen(Sh4 *sh4, Sh4OpArgs inst);

// DIV0S Rm, Rn
// 0010nnnnmmmm0111
void sh4_inst_binary_div0s_gen_gen(Sh4 *sh4, Sh4OpArgs inst);

// DIV0U
// 0000000000011001
void sh4_inst_noarg_div0u(Sh4 *sh4, Sh4OpArgs inst);

// DMULS.L Rm, Rn
//0011nnnnmmmm1101
void sh4_inst_binary_dmulsl_gen_gen(Sh4 *sh4, Sh4OpArgs inst);

// DMULU.L Rm, Rn
// 0011nnnnmmmm0101
void sh4_inst_binary_dmulul_gen_gen(Sh4 *sh4, Sh4OpArgs inst);

// EXTS.B Rm, Rn
// 0110nnnnmmmm1110
void sh4_inst_binary_extsb_gen_gen(Sh4 *sh4, Sh4OpArgs inst);

// EXTS.W Rm, Rnn
// 0110nnnnmmmm1111
void sh4_inst_binary_extsw_gen_gen(Sh4 *sh4, Sh4OpArgs inst);

// EXTU.B Rm, Rn
// 0110nnnnmmmm1100
void sh4_inst_binary_extub_gen_gen(Sh4 *sh4, Sh4OpArgs inst);

// EXTU.W Rm, Rn
// 0110nnnnmmmm1101
void sh4_inst_binary_extuw_gen_gen(Sh4 *sh4, Sh4OpArgs inst);

// MUL.L Rm, Rn
// 0000nnnnmmmm0111
void sh4_inst_binary_mull_gen_gen(Sh4 *sh4, Sh4OpArgs inst);

// MULS.W Rm, Rn
// 0010nnnnmmmm1111
void sh4_inst_binary_mulsw_gen_gen(Sh4 *sh4, Sh4OpArgs inst);

// MULU.W Rm, Rn
// 0010nnnnmmmm1110
void sh4_inst_binary_muluw_gen_gen(Sh4 *sh4, Sh4OpArgs inst);

// NEG Rm, Rn
// 0110nnnnmmmm1011
void sh4_inst_binary_neg_gen_gen(Sh4 *sh4, Sh4OpArgs inst);

// NEGC Rm, Rn
// 0110nnnnmmmm1010
void sh4_inst_binary_negc_gen_gen(Sh4 *sh4, Sh4OpArgs inst);

// SUB Rm, Rn
// 0011nnnnmmmm1000
void sh4_inst_binary_sub_gen_gen(Sh4 *sh4, Sh4OpArgs inst);

// SUBC Rm, Rn
// 0011nnnnmmmm1010
void sh4_inst_binary_subc_gen_gen(Sh4 *sh4, Sh4OpArgs inst);

// SUBV Rm, Rn
// 0011nnnnmmmm1011
void sh4_inst_binary_subv_gen_gen(Sh4 *sh4, Sh4OpArgs inst);

// AND Rm, Rn
// 0010nnnnmmmm1001
void sh4_inst_binary_and_gen_gen(Sh4 *sh4, Sh4OpArgs inst);

// NOT Rm, Rn
// 0110nnnnmmmm0111
void sh4_inst_binary_not_gen_gen(Sh4 *sh4, Sh4OpArgs inst);

// OR Rm, Rn
// 0010nnnnmmmm1011
void sh4_inst_binary_or_gen_gen(Sh4 *sh4, Sh4OpArgs inst);

// TST Rm, Rn
// 0010nnnnmmmm1000
void sh4_inst_binary_tst_gen_gen(Sh4 *sh4, Sh4OpArgs inst);

// XOR Rm, Rn
// 0010nnnnmmmm1010
void sh4_inst_binary_xor_gen_gen(Sh4 *sh4, Sh4OpArgs inst);

// SHAD Rm, Rn
// 0100nnnnmmmm1100
void sh4_inst_binary_shad_gen_gen(Sh4 *sh4, Sh4OpArgs inst);

// SHLD Rm, Rn
// 0100nnnnmmmm1101
void sh4_inst_binary_shld_gen_gen(Sh4 *sh4, Sh4OpArgs inst);

// LDC Rm, Rn_BANK
// 0100mmmm1nnn1110
void sh4_inst_binary_ldc_gen_bank(Sh4 *sh4, Sh4OpArgs inst);

// LDC.L @Rm+, Rn_BANK
// 0100mmmm1nnn0111
void sh4_inst_binary_ldcl_indgeninc_bank(Sh4 *sh4, Sh4OpArgs inst);

// STC Rm_BANK, Rn
// 0000nnnn1mmm0010
void sh4_inst_binary_stc_bank_gen(Sh4 *sh4, Sh4OpArgs inst);

// STC.L Rm_BANK, @-Rn
// 0100nnnn1mmm0011
void sh4_inst_binary_stcl_bank_inddecgen(Sh4 *sh4, Sh4OpArgs inst);

// LDS Rm,MACH
// 0100mmmm00001010
void sh4_inst_binary_lds_gen_mach(Sh4 *sh4, Sh4OpArgs inst);

// LDS Rm, MACL
// 0100mmmm00011010
void sh4_inst_binary_lds_gen_macl(Sh4 *sh4, Sh4OpArgs inst);

// STS MACH, Rn
// 0000nnnn00001010
void sh4_inst_binary_sts_mach_gen(Sh4 *sh4, Sh4OpArgs inst);

// STS MACL, Rn
// 0000nnnn00011010
void sh4_inst_binary_sts_macl_gen(Sh4 *sh4, Sh4OpArgs inst);

// LDS Rm, PR
// 0100mmmm00101010
void sh4_inst_binary_lds_gen_pr(Sh4 *sh4, Sh4OpArgs inst);

// STS PR, Rn
// 0000nnnn00101010
void sh4_inst_binary_sts_pr_gen(Sh4 *sh4, Sh4OpArgs inst);

// LDS.L @Rm+, MACH
// 0100mmmm00000110
void sh4_inst_binary_ldsl_indgeninc_mach(Sh4 *sh4, Sh4OpArgs inst);

// LDS.L @Rm+, MACL
// 0100mmmm00010110
void sh4_inst_binary_ldsl_indgeninc_macl(Sh4 *sh4, Sh4OpArgs inst);

// STS.L MACH, @-Rn
// 0100mmmm00000010
void sh4_inst_binary_stsl_mach_inddecgen(Sh4 *sh4, Sh4OpArgs inst);

// STS.L MACL, @-Rn
// 0100mmmm00010010
void sh4_inst_binary_stsl_macl_inddecgen(Sh4 *sh4, Sh4OpArgs inst);

// LDS.L @Rm+, PR
// 0100mmmm00100110
void sh4_inst_binary_ldsl_indgeninc_pr(Sh4 *sh4, Sh4OpArgs inst);

// STS.L PR, @-Rn
// 0100nnnn00100010
void sh4_inst_binary_stsl_pr_inddecgen(Sh4 *sh4, Sh4OpArgs inst);

// MOV.B Rm, @Rn
// 0010nnnnmmmm0000
void sh4_inst_binary_movb_gen_indgen(Sh4 *sh4, Sh4OpArgs inst);

// MOV.W Rm, @Rn
// 0010nnnnmmmm0001
void sh4_inst_binary_movw_gen_indgen(Sh4 *sh4, Sh4OpArgs inst);

// MOV.L Rm, @Rn
// 0010nnnnmmmm0010
void sh4_inst_binary_movl_gen_indgen(Sh4 *sh4, Sh4OpArgs inst);

// MOV.B @Rm, Rn
// 0110nnnnmmmm0000
void sh4_inst_binary_movb_indgen_gen(Sh4 *sh4, Sh4OpArgs inst);

// MOV.W @Rm, Rn
// 0110nnnnmmmm0001
void sh4_inst_binary_movw_indgen_gen(Sh4 *sh4, Sh4OpArgs inst);

// MOV.L @Rm, Rn
// 0110nnnnmmmm0010
void sh4_inst_binary_movl_indgen_gen(Sh4 *sh4, Sh4OpArgs inst);

// MOV.B Rm, @-Rn
// 0010nnnnmmmm0100
void sh4_inst_binary_movb_gen_inddecgen(Sh4 *sh4, Sh4OpArgs inst);

// MOV.W Rm, @-Rn
// 0010nnnnmmmm0101
void sh4_inst_binary_movw_gen_inddecgen(Sh4 *sh4, Sh4OpArgs inst);

// MOV.L Rm, @-Rn
// 0010nnnnmmmm0110
void sh4_inst_binary_movl_gen_inddecgen(Sh4 *sh4, Sh4OpArgs inst);

// MOV.B @Rm+, Rn
// 0110nnnnmmmm0100
void sh4_inst_binary_movb_indgeninc_gen(Sh4 *sh4, Sh4OpArgs inst);

// MOV.W @Rm+, Rn
// 0110nnnnmmmm0101
void sh4_inst_binary_movw_indgeninc_gen(Sh4 *sh4, Sh4OpArgs inst);

// MOV.L @Rm+, Rn
// 0110nnnnmmmm0110
void sh4_inst_binary_movl_indgeninc_gen(Sh4 *sh4, Sh4OpArgs inst);

// MAC.L @Rm+, @Rn+
// 0000nnnnmmmm1111
void sh4_inst_binary_macl_indgeninc_indgeninc(Sh4 *sh4, Sh4OpArgs inst);

// MAC.W @Rm+, @Rn+
// 0100nnnnmmmm1111
void sh4_inst_binary_macw_indgeninc_indgeninc(Sh4 *sh4, Sh4OpArgs inst);

// MOV.B R0, @(disp, Rn)
// 10000000nnnndddd
void sh4_inst_binary_movb_r0_binind_disp_gen(Sh4 *sh4, Sh4OpArgs inst);

// MOV.W R0, @(disp, Rn)
// 10000001nnnndddd
void sh4_inst_binary_movw_r0_binind_disp_gen(Sh4 *sh4, Sh4OpArgs inst);

// MOV.L Rm, @(disp, Rn)
// 0001nnnnmmmmdddd
void sh4_inst_binary_movl_gen_binind_disp_gen(Sh4 *sh4, Sh4OpArgs inst);

// MOV.B @(disp, Rm), R0
// 10000100mmmmdddd
void sh4_inst_binary_movb_binind_disp_gen_r0(Sh4 *sh4, Sh4OpArgs inst);

// MOV.W @(disp, Rm), R0
// 10000101mmmmdddd
void sh4_inst_binary_movw_binind_disp_gen_r0(Sh4 *sh4, Sh4OpArgs inst);

// MOV.L @(disp, Rm), Rn
// 0101nnnnmmmmdddd
void sh4_inst_binary_movl_binind_disp_gen_gen(Sh4 *sh4, Sh4OpArgs inst);

// MOV.B Rm, @(R0, Rn)
// 0000nnnnmmmm0100
void sh4_inst_binary_movb_gen_binind_r0_gen(Sh4 *sh4, Sh4OpArgs inst);

// MOV.W Rm, @(R0, Rn)
// 0000nnnnmmmm0101
void sh4_inst_binary_movw_gen_binind_r0_gen(Sh4 *sh4, Sh4OpArgs inst);

// MOV.L Rm, @(R0, Rn)
// 0000nnnnmmmm0110
void sh4_inst_binary_movl_gen_binind_r0_gen(Sh4 *sh4, Sh4OpArgs inst);

// MOV.B @(R0, Rm), Rn
// 0000nnnnmmmm1100
void sh4_inst_binary_movb_binind_r0_gen_gen(Sh4 *sh4, Sh4OpArgs inst);

// MOV.W @(R0, Rm), Rn
// 0000nnnnmmmm1101
void sh4_inst_binary_movw_binind_r0_gen_gen(Sh4 *sh4, Sh4OpArgs inst);

// MOV.L @(R0, Rm), Rn
// 0000nnnnmmmm1110
void sh4_inst_binary_movl_binind_r0_gen_gen(Sh4 *sh4, Sh4OpArgs inst);

// MOV.B R0, @(disp, GBR)
// 11000000dddddddd
void sh4_inst_binary_movb_r0_binind_disp_gbr(Sh4 *sh4, Sh4OpArgs inst);

// MOV.W R0, @(disp, GBR)
// 11000001dddddddd
void sh4_inst_binary_movw_r0_binind_disp_gbr(Sh4 *sh4, Sh4OpArgs inst);

// MOV.L R0, @(disp, GBR)
// 11000010dddddddd
void sh4_inst_binary_movl_r0_binind_disp_gbr(Sh4 *sh4, Sh4OpArgs inst);

// MOV.B @(disp, GBR), R0
// 11000100dddddddd
void sh4_inst_binary_movb_binind_disp_gbr_r0(Sh4 *sh4, Sh4OpArgs inst);

// MOV.W @(disp, GBR), R0
// 11000101dddddddd
void sh4_inst_binary_movw_binind_disp_gbr_r0(Sh4 *sh4, Sh4OpArgs inst);

// MOV.L @(disp, GBR), R0
// 11000110dddddddd
void sh4_inst_binary_movl_binind_disp_gbr_r0(Sh4 *sh4, Sh4OpArgs inst);

// MOVA @(disp, PC), R0
// 11000111dddddddd
void sh4_inst_binary_mova_binind_disp_pc_r0(Sh4 *sh4, Sh4OpArgs inst);

// MOVCA.L R0, @Rn
// 0000nnnn11000011
void sh4_inst_binary_movcal_r0_indgen(Sh4 *sh4, Sh4OpArgs inst);

// FLDI0 FRn
// 1111nnnn10001101
void sh4_inst_unary_fldi0_fr(Sh4 *sh4, Sh4OpArgs inst);

// FLDI1 Frn
// 1111nnnn10011101
void sh4_inst_unary_fldi1_fr(Sh4 *sh4, Sh4OpArgs inst);

// FMOV FRm, FRn
// 1111nnnnmmmm1100
void sh4_inst_binary_fmov_fr_fr(Sh4 *sh4, Sh4OpArgs inst);

// FMOV.S @Rm, FRn
// 1111nnnnmmmm1000
void sh4_inst_binary_fmovs_indgen_fr(Sh4 *sh4, Sh4OpArgs inst);

// FMOV.S @(R0,Rm), FRn
// 1111nnnnmmmm0110
void sh4_inst_binary_fmovs_binind_r0_gen_fr(Sh4 *sh4, Sh4OpArgs inst);

// FMOV.S @Rm+, FRn
// 1111nnnnmmmm1001
void sh4_inst_binary_fmovs_indgeninc_fr(Sh4 *sh4, Sh4OpArgs inst);

// FMOV.S FRm, @Rn
// 1111nnnnmmmm1010
void sh4_inst_binary_fmovs_fr_indgen(Sh4 *sh4, Sh4OpArgs inst);

// FMOV.S FRm, @-Rn
// 1111nnnnmmmm1011
void sh4_inst_binary_fmovs_fr_inddecgen(Sh4 *sh4, Sh4OpArgs inst);

// FMOV.S FRm, @(R0, Rn)
// 1111nnnnmmmm0111
void sh4_inst_binary_fmovs_fr_binind_r0_gen(Sh4 *sh4, Sh4OpArgs inst);

// FMOV DRm, DRn
// 1111nnn0mmm01100
void sh4_inst_binary_fmov_dr_dr(Sh4 *sh4, Sh4OpArgs inst);

// FMOV @Rm, DRn
// 1111nnn0mmmm1000
void sh4_inst_binary_fmov_indgen_dr(Sh4 *sh4, Sh4OpArgs inst);

// FMOV @(R0, Rm), DRn
// 1111nnn0mmmm0110
void sh4_inst_binary_fmov_binind_r0_gen_dr(Sh4 *sh4, Sh4OpArgs inst);

// FMOV @Rm+, DRn
// 1111nnn0mmmm1001
void sh4_inst_binary_fmov_indgeninc_dr(Sh4 *sh4, Sh4OpArgs inst);

// FMOV DRm, @Rn
// 1111nnnnmmm01010
void sh4_inst_binary_fmov_dr_indgen(Sh4 *sh4, Sh4OpArgs inst);

// FMOV DRm, @-Rn
// 1111nnnnmmm01011
void sh4_inst_binary_fmov_dr_inddecgen(Sh4 *sh4, Sh4OpArgs inst);

// FMOV DRm, @(R0,Rn)
// 1111nnnnmmm00111
void sh4_inst_binary_fmov_dr_binind_r0_gen(Sh4 *sh4, Sh4OpArgs inst);

// FLDS FRm, FPUL
// 1111mmmm00011101
void sh4_inst_binary_flds_fr_fpul(Sh4 *sh4, Sh4OpArgs inst);

// FSTS FPUL, FRn
// 1111nnnn00001101
void sh4_inst_binary_fsts_fpul_fr(Sh4 *sh4, Sh4OpArgs inst);

// FABS FRn
// 1111nnnn01011101
void sh4_inst_unary_fabs_fr(Sh4 *sh4, Sh4OpArgs inst);

// FADD FRm, FRn
// 1111nnnnmmmm0000
void sh4_inst_binary_fadd_fr_fr(Sh4 *sh4, Sh4OpArgs inst);

// FCMP/EQ FRm, FRn
// 1111nnnnmmmm0100
void sh4_inst_binary_fcmpeq_fr_fr(Sh4 *sh4, Sh4OpArgs inst);

// FCMP/GT FRm, FRn
// 1111nnnnmmmm0101
void sh4_inst_binary_fcmpgt_fr_fr(Sh4 *sh4, Sh4OpArgs inst);

// FDIV FRm, FRn
// 1111nnnnmmmm0011
void sh4_inst_binary_fdiv_fr_fr(Sh4 *sh4, Sh4OpArgs inst);

// FLOAT FPUL, FRn
// 1111nnnn00101101
void sh4_inst_binary_float_fpul_fr(Sh4 *sh4, Sh4OpArgs inst);

// FMAC FR0, FRm, FRn
// 1111nnnnmmmm1110
void sh4_inst_trinary_fmac_fr0_fr_fr(Sh4 *sh4, Sh4OpArgs inst);

// FMUL FRm, FRn
// 1111nnnnmmmm0010
void sh4_inst_binary_fmul_fr_fr(Sh4 *sh4, Sh4OpArgs inst);

// FNEG FRn
// 1111nnnn01001101
void sh4_inst_unary_fneg_fr(Sh4 *sh4, Sh4OpArgs inst);

// FSQRT FRn
// 1111nnnn01101101
void sh4_inst_unary_fsqrt_fr(Sh4 *sh4, Sh4OpArgs inst);

// FSUB FRm, FRn
// 1111nnnnmmmm0001
void sh4_inst_binary_fsub_fr_fr(Sh4 *sh4, Sh4OpArgs inst);

// FTRC FRm, FPUL
// 1111mmmm00111101
void sh4_inst_binary_ftrc_fr_fpul(Sh4 *sh4, Sh4OpArgs inst);

// FABS DRn
// 1111nnn001011101
void sh4_inst_unary_fabs_dr(Sh4 *sh4, Sh4OpArgs inst);

// FADD DRm, DRn
// 1111nnn0mmm00000
void sh4_inst_binary_fadd_dr_dr(Sh4 *sh4, Sh4OpArgs inst);

// FCMP/EQ DRm, DRn
// 1111nnn0mmm00100
void sh4_inst_binary_fcmpeq_dr_dr(Sh4 *sh4, Sh4OpArgs inst);

// FCMP/GT DRm, DRn
// 1111nnn0mmm00101
void sh4_inst_binary_fcmpgt_dr_dr(Sh4 *sh4, Sh4OpArgs inst);

// FDIV DRm, DRn
// 1111nnn0mmm00011
void sh4_inst_binary_fdiv_dr_dr(Sh4 *sh4, Sh4OpArgs inst);

// FCNVDS DRm, FPUL
// 1111mmm010111101
void sh4_inst_binary_fcnvds_dr_fpul(Sh4 *sh4, Sh4OpArgs inst);

// FCNVSD FPUL, DRn
// 1111nnn010101101
void sh4_inst_binary_fcnvsd_fpul_dr(Sh4 *sh4, Sh4OpArgs inst);

// FLOAT FPUL, DRn
// 1111nnn000101101
void sh4_inst_binary_float_fpul_dr(Sh4 *sh4, Sh4OpArgs inst);

// FMUL DRm, DRn
// 1111nnn0mmm00010
void sh4_inst_binary_fmul_dr_dr(Sh4 *sh4, Sh4OpArgs inst);

// FNEG DRn
// 1111nnn001001101
void sh4_inst_unary_fneg_dr(Sh4 *sh4, Sh4OpArgs inst);

// FSQRT DRn
// 1111nnn001101101
void sh4_inst_unary_fsqrt_dr(Sh4 *sh4, Sh4OpArgs inst);

// FSUB DRm, DRn
// 1111nnn0mmm00001
void sh4_inst_binary_fsub_dr_dr(Sh4 *sh4, Sh4OpArgs inst);

// FTRC DRm, FPUL
// 1111mmm000111101
void sh4_inst_binary_ftrc_dr_fpul(Sh4 *sh4, Sh4OpArgs inst);

// LDS Rm, FPSCR
// 0100mmmm01101010
void sh4_inst_binary_lds_gen_fpscr(Sh4 *sh4, Sh4OpArgs inst);

// LDS Rm, FPUL
// 0100mmmm01011010
void sh4_inst_binary_gen_fpul(Sh4 *sh4, Sh4OpArgs inst);

// LDS.L @Rm+, FPSCR
// 0100mmmm01100110
void sh4_inst_binary_ldsl_indgeninc_fpscr(Sh4 *sh4, Sh4OpArgs inst);

// LDS.L @Rm+, FPUL
// 0100mmmm01010110
void sh4_inst_binary_ldsl_indgeninc_fpul(Sh4 *sh4, Sh4OpArgs inst);

// STS FPSCR, Rn
// 0000nnnn01101010
void sh4_inst_binary_sts_fpscr_gen(Sh4 *sh4, Sh4OpArgs inst);

// STS FPUL, Rn
// 0000nnnn01011010
void sh4_inst_binary_sts_fpul_gen(Sh4 *sh4, Sh4OpArgs inst);

// STS.L FPSCR, @-Rn
// 0100nnnn01100010
void sh4_inst_binary_stsl_fpscr_inddecgen(Sh4 *sh4, Sh4OpArgs inst);

// STS.L FPUL, @-Rn
// 0100nnnn01010010
void sh4_inst_binary_stsl_fpul_inddecgen(Sh4 *sh4, Sh4OpArgs inst);

// FMOV DRm, XDn
// 1111nnn1mmm01100
void sh4_inst_binary_fmov_dr_xd(Sh4 *sh4, Sh4OpArgs inst);

// FMOV XDm, DRn
// 1111nnn0mmm11100
void sh4_inst_binary_fmov_xd_dr(Sh4 *sh4, Sh4OpArgs inst);

// FMOV XDm, XDn
// 1111nnn1mmm11100
void sh4_inst_binary_fmov_xd_xd(Sh4 *sh4, Sh4OpArgs inst);

// FMOV @Rm, XDn
// 1111nnn1mmmm1000
void sh4_inst_binary_fmov_indgen_xd(Sh4 *sh4, Sh4OpArgs inst);

// FMOV @Rm+, XDn
// 1111nnn1mmmm1001
void sh4_inst_binary_fmov_indgeninc_xd(Sh4 *sh4, Sh4OpArgs inst);

// FMOV @(R0, Rn), XDn
// 1111nnn1mmmm0110
void sh4_inst_binary_fmov_binind_r0_gen_xd(Sh4 *sh4, Sh4OpArgs inst);

// FMOV XDm, @Rn
// 1111nnnnmmm11010
void sh4_inst_binary_fmov_xd_indgen(Sh4 *sh4, Sh4OpArgs inst);

// FMOV XDm, @-Rn
// 1111nnnnmmm11011
void sh4_inst_binary_fmov_xd_inddecgen(Sh4 *sh4, Sh4OpArgs inst);

// FMOV XDm, @(R0, Rn)
// 1111nnnnmmm10111
void sh4_inst_binary_fmov_xs_binind_r0_gen(Sh4 *sh4, Sh4OpArgs inst);

// FIPR FVm, FVn - vector dot product
// 1111nnmm11101101
void sh4_inst_binary_fipr_fv_fv(Sh4 *sh4, Sh4OpArgs inst);

// FTRV MXTRX, FVn - multiple vector by matrix
// 1111nn0111111101
void sh4_inst_binary_fitrv_mxtrx_fv(Sh4 *sh4, Sh4OpArgs inst);

// FSRRA FRn
// 1111nnnn01111101
void sh4_inst_unary_fsrra_frn(Sh4  *sh4, Sh4OpArgs inst);

/*
 * fake opcode used for implementing softbreaks and errors.
 *
 * Depending on preprocessor options, this function might signal to the
 * debugger that there was a TRAP, or it might throw a C++ function or it
 * might set an sh4 CPU exception (which is what would happen on an actual
 * sh4).
 */
void sh4_inst_invalid(Sh4 *sh4, Sh4OpArgs inst);

////////////////////////////////////////////////////////////////////////////////
//
// The following handlers are for floating-point opcodes that share their
// opcodes with other floating-point opcodes.  which handler gets called is
// based on either the PR bit or the SZ bit in the FPSCR
//
////////////////////////////////////////////////////////////////////////////////

#define FPU_HANDLER(name) sh4_fpu_inst_ ## name

#define DECL_FPU_HANDLER(name) void FPU_HANDLER(name) (Sh4 *sh4, Sh4OpArgs inst)

#define DEF_FPU_HANDLER(name, mask, on_false, on_true) \
    DECL_FPU_HANDLER(name) {                           \
        if (sh4->reg[SH4_REG_FPSCR] & (mask)) {        \
            on_true(sh4, inst);                        \
        } else {                                       \
            on_false(sh4, inst);                       \
        }                                              \
    }

#define DEF_FPU_HANDLER_CUSTOM(name)            \
    DECL_FPU_HANDLER(name)

// these handlers have no specified behavior for when the PR bit is set

// FLDI0 FRn
// 1111nnnn10001101
DECL_FPU_HANDLER(fldi0);

// FLDI1 Frn
// 1111nnnn10011101
DECL_FPU_HANDLER(fldi1);

// these handlers depend on the SZ bit
// FMOV FRm, FRn
// 1111nnnnmmmm1100
// FMOV DRm, DRn
// 1111nnn0mmm01100
DECL_FPU_HANDLER(fmov_gen);

// FMOV.S @Rm, FRn
// 1111nnnnmmmm1000
// FMOV @Rm, DRn
// 1111nnn0mmmm1000
DECL_FPU_HANDLER(fmovs_ind_gen);

// FMOV.S @(R0, Rm), FRn
// 1111nnnnmmmm0110
// FMOV @(R0, Rm), DRn
// 1111nnn0mmmm0110
DECL_FPU_HANDLER(fmov_binind_r0_gen_fpu);

// FMOV.S @Rm+, FRn
// 1111nnnnmmmm1001
// FMOV @Rm+, DRn
// 1111nnn0mmmm1001
DECL_FPU_HANDLER(fmov_indgeninc_fpu);

// FMOV.S FRm, @Rn
// 1111nnnnmmmm1010
// FMOV DRm, @Rn
// 1111nnnnmmm01010
DECL_FPU_HANDLER(fmov_fpu_indgen);

// FMOV.S FRm, @-Rn
// 1111nnnnmmmm1011
// FMOV DRm, @-Rn
// 1111nnnnmmm01011
DECL_FPU_HANDLER(fmov_fpu_inddecgen);

// FMOV.S FRm, @(R0, Rn)
// 1111nnnnmmmm0111
// FMOV DRm, @(R0, Rn)
// 1111nnnnmmm00111
DECL_FPU_HANDLER(fmov_fpu_binind_r0_gen);

// FABS FRn
// 1111nnnn01011101
// FABS DRn
// 1111nnn001011101
DECL_FPU_HANDLER(fabs_fpu);

// FADD FRm, FRn
// 1111nnnnmmmm0000
// FADD DRm, DRn
// 1111nnn0mmm00000
DECL_FPU_HANDLER(fadd_fpu);

// FCMP/EQ FRm, FRn
// 1111nnnnmmmm0100
// FCMP/EQ DRm, DRn
// 1111nnn0mmm00100
DECL_FPU_HANDLER(fcmpeq_fpu);

// FCMP/GT FRm, FRn
// 1111nnnnmmmm0101
// FCMP/GT DRm, DRn
// 1111nnn0mmm00101
DECL_FPU_HANDLER(fcmpgt_fpu);

// FDIV FRm, FRn
// 1111nnnnmmmm0011
// FDIV DRm, DRn
// 1111nnn0mmm00011
DECL_FPU_HANDLER(fdiv_fpu);

// FLOAT FPUL, FRn
// 1111nnnn00101101
// FLOAT FPUL, DRn
// 1111nnn000101101
DECL_FPU_HANDLER(float_fpu);

// FMAC FR0, FRm, FRn
// 1111nnnnmmmm1110
DECL_FPU_HANDLER(fmac_fpu);

// FMUL FRm, FRn
// 1111nnnnmmmm0010
// FMUL DRm, DRn
// 1111nnn0mmm00010
DECL_FPU_HANDLER(fmul_fpu);

// FNEG FRn
// 1111nnnn01001101
// FNEG DRn
// 1111nnn001001101
DECL_FPU_HANDLER(fneg_fpu);

// FSQRT FRn
// 1111nnnn01101101
// FSQRT DRn
// 1111nnn001101101
DECL_FPU_HANDLER(fsqrt_fpu);

// FSUB FRm, FRn
// 1111nnnnmmmm0001
// FSUB DRm, DRn
// 1111nnn0mmm00001
DECL_FPU_HANDLER(fsub_fpu);

// FTRC FRm, FPUL
// 1111mmmm00111101
// FTRC DRm, FPUL
// 1111mmm000111101
DECL_FPU_HANDLER(ftrc_fpu);

// FCNVDS DRm, FPUL
// 1111mmm010111101
DECL_FPU_HANDLER(fcnvds_fpu);

// FCNVSD FPUL, DRn
// 1111nnn010101101
DECL_FPU_HANDLER(fcnvsd_fpu);

// FSCA FPUL, DRn
// 1111nnn011111101
DECL_FPU_HANDLER(fsca_fpu);

// FSRRA FRn
// 1111nnnn01111101
DECL_FPU_HANDLER(fsrra_fpu);

#endif
