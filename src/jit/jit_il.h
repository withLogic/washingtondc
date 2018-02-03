/*******************************************************************************
 *
 *
 *    WashingtonDC Dreamcast Emulator
 *    Copyright (C) 2018 snickerbockers
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

#ifndef JIT_IL_H_
#define JIT_IL_H_

// for union Sh4OpArgs
#include "hw/sh4/sh4_inst.h"

enum jit_opcode {
    // this opcode calls an interpreter function
    JIT_OP_FALLBACK,

    /*
     * this stores a given register plus a constant
     * offset as the branch destination
     */
    JIT_OP_PREPARE_JUMP,

    // this sets a given constant as the branch destination
    JIT_OP_PREPARE_JUMP_CONST,

    /*
     * This stores a given register plus a constant offset as a jump
     * destination for a failed conditional jump.
     */
    JIT_OP_PREPARE_ALT_JUMP,

    // This jumps to the jump destination address previously stored
    JIT_OP_JUMP,

    /*
     * This can be configured to set the conditional jump flag if t is set, or
     * it can be configured to set the conditional jump flag if t is not set.
     */
    JIT_SET_COND_JUMP_BASED_ON_T,

    // this will jump iff the conditional jump flag is set
    JIT_JUMP_COND,

    // this will set a register to the given constant value
    JIT_SET_SLOT,

    // this will copy a slot into SR and handle any state changes
    JIT_OP_RESTORE_SR,

    // read 16 bits from a constant address and store them in a given slot
    JIT_OP_READ_16_SLOT,

    // sign-extend a 16-bit int in a register into a 32-bit int
    JIT_OP_SIGN_EXTEND_16,

    JIT_OP_READ_32_SLOT,

    // load 32-bits from a host memory address into a jit register
    JIT_OP_LOAD_SLOT,

    // store 32-bits from a jit register address into a host memory address
    JIT_OP_STORE_SLOT,

    // add one slot into another
    JIT_OP_ADD
};

struct jit_fallback_immed {
    void(*fallback_fn)(Sh4*,Sh4OpArgs);
    Sh4OpArgs inst;
};

struct prepare_jump_immed {
    unsigned slot_idx;
    unsigned offs; // constant offset added to the register
};

struct prepare_jump_const_immed {
    unsigned new_pc;
};

struct prepare_alt_jump_immed {
    unsigned new_pc;
};

struct set_cond_jump_based_on_t_immed {
    unsigned t_flag;
};

struct set_slot_immed {
    unsigned slot_idx;
    uint32_t new_val;
};

struct restore_sr_immed {
    unsigned slot_no;
};

struct read_16_slot_immed {
    addr32_t addr;
    unsigned slot_no;
};

struct sign_extend_16_immed {
    unsigned slot_no;
};

struct read_32_slot_immed {
    addr32_t addr;
    unsigned slot_no;
};

struct load_slot_immed {
    uint32_t const *src;
    unsigned slot_no;
};

struct store_slot_immed {
    uint32_t *dst;
    unsigned slot_no;
};

struct add_immed {
    unsigned slot_src, slot_dst;
};

union jit_immed {
    struct jit_fallback_immed fallback;
    struct prepare_jump_immed prepare_jump;
    struct prepare_jump_const_immed prepare_jump_const;
    struct prepare_alt_jump_immed prepare_alt_jump;
    struct set_cond_jump_based_on_t_immed set_cond_jump_based_on_t;
    struct set_slot_immed set_slot;
    struct restore_sr_immed restore_sr;
    struct read_16_slot_immed read_16_slot;
    struct sign_extend_16_immed sign_extend_16;
    struct read_32_slot_immed read_32_slot;
    struct load_slot_immed load_slot;
    struct store_slot_immed store_slot;
    struct add_immed add;
};

struct jit_inst {
    enum jit_opcode op;
    union jit_immed immed;
};

void jit_fallback(struct jit_inst *op,
                  void(*fallback_fn)(Sh4*,Sh4OpArgs), inst_t inst);
void jit_prepare_jump(struct jit_inst *op, unsigned slot_idx, unsigned offs);
void jit_prepare_jump_const(struct jit_inst *op, unsigned new_pc);
void jit_prepare_alt_jump(struct jit_inst *op, unsigned new_pc);
void jit_jump(struct jit_inst *op);
void jit_set_cond_jump_based_on_t(struct jit_inst *op, unsigned t_val);
void jit_jump_cond(struct jit_inst *op);
void jit_set_slot(struct jit_inst *op, unsigned reg_idx, uint32_t new_val);
void jit_restore_sr(struct jit_inst *op, unsigned slot_no);
void jit_read_16_slot(struct jit_inst *op, addr32_t addr, unsigned slot_no);
void jit_sign_extend_16(struct jit_inst *op, unsigned slot_no);
void jit_read_32_slot(struct jit_inst *op, addr32_t addr, unsigned slot_no);
void jit_load_slot(struct jit_inst *op, unsigned slot_no, uint32_t const *src);
void jit_store_slot(struct jit_inst *op, unsigned slot_no, uint32_t *dst);
void jit_add(struct jit_inst *op, unsigned slot_src, unsigned slot_dst);

#endif
