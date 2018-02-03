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

#include "jit_il.h"

void jit_fallback(struct jit_inst *op,
                  void(*fallback_fn)(Sh4*,Sh4OpArgs), inst_t inst) {
    op->op = JIT_OP_FALLBACK;
    op->immed.fallback.fallback_fn = fallback_fn;
    op->immed.fallback.inst.inst = inst;
}

void jit_prepare_jump(struct jit_inst *op, unsigned slot_idx) {
    op->op = JIT_OP_PREPARE_JUMP;
    op->immed.prepare_jump.slot_idx = slot_idx;
}

void jit_prepare_jump_const(struct jit_inst *op, unsigned new_pc) {
    op->op = JIT_OP_PREPARE_JUMP_CONST;
    op->immed.prepare_jump_const.new_pc = new_pc;
}

void jit_prepare_alt_jump(struct jit_inst *op, unsigned new_pc) {
    op->op = JIT_OP_PREPARE_ALT_JUMP;
    op->immed.prepare_alt_jump.new_pc = new_pc;
}

void jit_jump(struct jit_inst *op) {
    op->op = JIT_OP_JUMP;
}

void jit_set_cond_jump_based_on_t(struct jit_inst *op, unsigned t_val) {
    op->op = JIT_SET_COND_JUMP_BASED_ON_T;
    op->immed.set_cond_jump_based_on_t.t_flag = t_val;
}

void jit_jump_cond(struct jit_inst *op) {
    op->op = JIT_JUMP_COND;
}

void jit_set_slot(struct jit_inst *op, unsigned slot_idx, uint32_t new_val) {
    op->op = JIT_SET_SLOT;
    op->immed.set_slot.new_val = new_val;
    op->immed.set_slot.slot_idx = slot_idx;
}

void jit_restore_sr(struct jit_inst *op, unsigned slot_no) {
    op->op = JIT_OP_RESTORE_SR;
    op->immed.restore_sr.slot_no = slot_no;
}

void jit_read_16_slot(struct jit_inst *op, addr32_t addr, unsigned slot_no) {
    op->op = JIT_OP_READ_16_SLOT;
    op->immed.read_16_slot.addr = addr;
    op->immed.read_16_slot.slot_no = slot_no;
}

void jit_sign_extend_16(struct jit_inst *op, unsigned slot_no) {
    op->op = JIT_OP_SIGN_EXTEND_16;
    op->immed.sign_extend_16.slot_no = slot_no;
}

void jit_read_32_slot(struct jit_inst *op, addr32_t addr, unsigned slot_no) {
    op->op = JIT_OP_READ_32_SLOT;
    op->immed.read_32_slot.addr = addr;
    op->immed.read_32_slot.slot_no = slot_no;
}

void jit_load_slot(struct jit_inst *op, unsigned slot_no, uint32_t const *src) {
    op->op = JIT_OP_LOAD_SLOT;
    op->immed.load_slot.src = src;
    op->immed.load_slot.slot_no = slot_no;
}

void jit_store_slot(struct jit_inst *op, unsigned slot_no, uint32_t *dst) {
    op->op = JIT_OP_STORE_SLOT;
    op->immed.store_slot.dst = dst;
    op->immed.store_slot.slot_no = slot_no;
}

void jit_add(struct jit_inst *op, unsigned slot_src, unsigned slot_dst) {
    op->op = JIT_OP_ADD;
    op->immed.add.slot_src = slot_src;
    op->immed.add.slot_dst = slot_dst;
}

void jit_add_const32(struct jit_inst *op, unsigned slot_dst, uint32_t const32) {
    op->op = JIT_OP_ADD_CONST32;
    op->immed.add_const32.slot_dst = slot_dst;
    op->immed.add_const32.const32 = const32;
}
