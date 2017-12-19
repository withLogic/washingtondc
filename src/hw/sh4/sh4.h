/*******************************************************************************
 *
 *
 *    WashingtonDC Dreamcast Emulator
 *    Copyright (C) 2016, 2017 snickerbockers
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

#ifndef SH4_H_
#define SH4_H_

#include <assert.h>
#include <stdint.h>
#include <stdatomic.h>

#ifdef ENABLE_SH4_MMU
#include "sh4_mmu.h"
#endif

#include "error.h"
#include "types.h"
#include "sh4_inst.h"
#include "sh4_reg.h"
#include "sh4_mem.h"
#include "sh4_tmu.h"
#include "sh4_ocache.h"
#include "sh4_excp.h"
#include "sh4_scif.h"
#include "sh4_dmac.h"

/*
 * The clock-scale is here defined as the number of scheduler cyclers per sh4
 * cycle.
 *
 * To convert dc_sched cycles to sh4 cycles, divide by SH4_CLOCK_SCALE
 * To convert sh4 cycles to dc_sched_cycles, multiply by SH4_CLOCK_SCALE
 */
#define SH4_CLOCK_SCALE (SCHED_FREQUENCY / (200 * 1000 * 1000))

static_assert(SCHED_FREQUENCY % (200 * 1000 * 1000) == 0,
              "scheduler frequency does not cleanly divide by SH4 frequency");

/* Hitachi SuperH-4 interpreter */

#define SH4_N_FLOAT_REGS 16
#define SH4_N_DOUBLE_REGS 8

enum Sh4ExecState {
    SH4_EXEC_STATE_NORM,
    SH4_EXEC_STATE_SLEEP,
    SH4_EXEC_STATE_STANDBY
};
typedef enum Sh4ExecState Sh4ExecState;

struct Sh4 {
    Sh4ExecState exec_state;

    reg32_t reg[SH4_REGISTER_COUNT];

#ifdef ENABLE_SH4_MMU
    struct sh4_mmu mmu;
#endif
    /*
     * If the CPU is executing a delayed branch instruction, then
     * delayed_branch will be true and delayed_branch_addr will point to the
     * address to branch to.  After executing one instruction, delayed_branch
     * will be set to false and the CPU will jump to delayed_branch_addr.
     *
     * If the branch instruction evaluates to false (ie, there is not a delayed
     * branch) then delayed_branch will never be set to true.  This means that
     * the interpreter will not raise any exceptions caused by executing a
     * branch instruction in a delay slot; this is an inaccuracy which may need
     * to be revisited in the future.
     */
    bool delayed_branch;
    addr32_t delayed_branch_addr;

    struct sh4_tmu tmu;

    /*
     * operand cache - this is really only here to be used as RAM
     * when the ORA bit is set in CCR
     */
    struct sh4_ocache ocache;

    struct sh4_intc intc;

    struct sh4_scif scif;

    struct sh4_dmac dmac;

    /*
     * pointer to place where memory-mapped registers are stored.
     * RegReadHandlers and RegWriteHandlers do not need to use this as long as
     * they are consistent.
     */
    uint8_t *reg_area;

    /*
     * this is used by sh4_count_inst_cycles to track the type of the last
     * instruction that was executed.  This is used to determine if the next
     * instruction to be executed should advance the cycle count, or if it
     * would have been executed by the second pipeline on a real sh4.
     */
     sh4_inst_group_t last_inst_type;

#ifdef ENABLE_DEBUGGER
    /*
     * this member is used to implement watchpoints.  When a watchpoint
     * is hit by sh4_write_mem or sh4_read_mem, this will be set to true
     * so that lower layers in the call-stack know the operation was aborted.
     * This is needed to handle watchpoints that happen in delayed-branch slots.
     */
    bool aborted_operation;
#endif
};

typedef struct Sh4 Sh4;

void sh4_init(Sh4 *sh4);
void sh4_cleanup(Sh4 *sh4);

// reset all values to their power-on-reset values
void sh4_on_hard_reset(Sh4 *sh4);

// returns the program counter
reg32_t sh4_get_pc(Sh4 *sh4);

/*
 * call this function instead of setting the value directly to make sure
 * that any state changes are immediately processed.
 */
void sh4_set_fpscr(Sh4 *sh4, reg32_t new_val);

// clear the cause bits in the FPSCR reg
static inline void sh4_fpu_clear_cause(Sh4 *sh4) {
#ifndef SH4_FPU_FAST
    sh4->reg[SH4_REG_FPSCR] &= ~SH4_FPSCR_CAUSE_MASK;
#endif
}

// these four APIs are intended primarily for debuggers to use
void sh4_get_regs(Sh4 *sh4, reg32_t reg_out[SH4_REGISTER_COUNT]);
/* FpuReg sh4_get_fpu(Sh4 *sh4); */
void sh4_set_regs(Sh4 *sh4, reg32_t const reg_out[SH4_REGISTER_COUNT]);
void sh4_set_individual_reg(Sh4 *sh4, unsigned reg_no, reg32_t reg_val);
/* void sh4_set_fpu(Sh4 *sh4, FpuReg src); */

void sh4_bank_switch(Sh4 *sh4);
void sh4_bank_switch_maybe(Sh4 *sh4, reg32_t old_sr, reg32_t new_sr);

void sh4_fpu_bank_switch(Sh4 *sh4);
void sh4_fpu_bank_switch_maybe(Sh4 *sh4, reg32_t old_fpscr, reg32_t new_fpscr);

static inline void sh4_next_inst(Sh4 *sh4) {
    sh4->reg[SH4_REG_PC] += 2;
}

/*
 * return the index of the given general-purpose register.
 * This function takes bank-switching into account.
 */
static inline sh4_reg_idx_t sh4_gen_reg_idx(Sh4 *sh4, int reg_no) {
    assert(!(reg_no & ~0xf));

    return (sh4_reg_idx_t)(SH4_REG_R0 + reg_no);
}

/*
 * return a pointer to the given general-purpose register.
 * This function takes bank-switching into account.
 */
static inline reg32_t *sh4_gen_reg(Sh4 *sh4, int idx) {
    return sh4->reg + sh4_gen_reg_idx(sh4, idx);
}

/* return an index to the given banked general-purpose register */
static inline sh4_reg_idx_t sh4_bank_reg_idx(Sh4 *sh4, int idx) {
    assert(!(idx & ~0x7));

    return (sh4_reg_idx_t)(SH4_REG_R0_BANK + idx);
}

// return a pointer to the given banked general-purpose register
static inline reg32_t *sh4_bank_reg(Sh4 *sh4, int idx) {
    return sh4->reg + sh4_bank_reg_idx(sh4, idx);
}

static inline sh4_reg_idx_t sh4_bank0_reg_idx(Sh4 *sh4, int idx) {
    assert(!(idx & ~0x7));

    if (sh4->reg[SH4_REG_SR] & SH4_SR_RB_MASK)
        return (sh4_reg_idx_t)(SH4_REG_R0_BANK + idx);
    return (sh4_reg_idx_t)(SH4_REG_R0 + idx);
}

static inline reg32_t *sh4_bank0_reg(Sh4 *sh4, int idx) {
    return sh4->reg + sh4_bank0_reg_idx(sh4, idx);
}

static inline sh4_reg_idx_t sh4_bank1_reg_idx(Sh4 *sh4, int idx) {
    assert(!(idx & ~0x7));

    if (sh4->reg[SH4_REG_SR] & SH4_SR_RB_MASK)
        return (sh4_reg_idx_t)(SH4_REG_R0 + idx);
    return (sh4_reg_idx_t)(SH4_REG_R0_BANK + idx);
}

static inline reg32_t *sh4_bank1_reg(Sh4 *sh4, int idx) {
    return sh4->reg + sh4_bank1_reg_idx(sh4, idx);
}

/*
 * access single-precision floating-point register,
 * taking bank-switching into account
 */
static inline float *sh4_fpu_fr(Sh4 *sh4, unsigned reg_no) {
    assert(reg_no < SH4_N_FLOAT_REGS);

    return (float*)(sh4->reg + SH4_REG_FR0 + reg_no);
}

static inline float *sh4_fpu_xf(Sh4 *sh4, unsigned reg_no) {
    assert(reg_no < SH4_N_FLOAT_REGS);

    return (float*)(sh4->reg + SH4_REG_XF0 + reg_no);
}

static inline float *sh4_bank0_fpu_fr(Sh4 *sh4, unsigned reg_no) {
    if (sh4->reg[SH4_REG_FPSCR] & SH4_FPSCR_FR_MASK)
        return sh4_fpu_xf(sh4, reg_no);
    return sh4_fpu_fr(sh4, reg_no);
}

static inline float *sh4_bank1_fpu_fr(Sh4 *sh4, unsigned reg_no) {
    if (sh4->reg[SH4_REG_FPSCR] & SH4_FPSCR_FR_MASK)
        return sh4_fpu_fr(sh4, reg_no);
    return sh4_fpu_xf(sh4, reg_no);
}

/*
 * access double-precision floating-point register,
 * taking bank-switching into account
 */
static inline double *sh4_fpu_dr(Sh4 *sh4, unsigned reg_no) {
    assert(reg_no < SH4_N_DOUBLE_REGS);

    return (double*)(sh4->reg + SH4_REG_FR0 + (reg_no << 1));
}

static inline double *sh4_fpu_xd(Sh4 *sh4, unsigned reg_no) {
    assert(reg_no < SH4_N_DOUBLE_REGS);

    return (double*)(sh4->reg + SH4_REG_XD0 + (reg_no << 1));
}

/*
 * this function should be called every time sr has just been written to and
 * bits other than T/Q/M/S may have changed
 */
void sh4_on_sr_change(Sh4 *sh4, reg32_t old_sr);

/*
 * The purpose of this function is to do things which need to be performed
 * periodically, but not with any urgency or hard timing requirements.
 *
 * Currently, that means the only thing it does is check to see if the serial
 * server wants to communicate with the SCIF; in the future other tasks may be
 * added in here as well if I need them.
 */
static inline void sh4_periodic(Sh4 *sh4) {
    if (!atomic_flag_test_and_set(&sh4->scif.nothing_pending))
        sh4_scif_periodic(sh4);
}

dc_cycle_stamp_t sh4_get_cycles(void);

// Fetches the given instruction's metadata and returns it.
static inline opcode_func_t
sh4_decode_inst(Sh4 *sh4, inst_t inst) {
    return sh4_inst_lut[inst];
}

/*
 * return the number of cycles this instruction requires.  This is not the same
 * as the instruction's issue cycles due to the dual-issue pipeline of the sh4.
 */
static inline unsigned
sh4_count_inst_cycles(Sh4 *sh4, sh4_inst_group_t group, unsigned issue) {
    unsigned n_cycles;
    if ((sh4->last_inst_type == SH4_GROUP_NONE) ||
        ((group == SH4_GROUP_CO) ||
         (sh4->last_inst_type == SH4_GROUP_CO) ||
         ((sh4->last_inst_type == group) && (group != SH4_GROUP_MT)))) {
        // This instruction was not free
        n_cycles = issue;

        /*
         * no need to check for SH4_GROUP_CO here because we'll do that when we
         * check for last_inst_type==SH4_GROUP_CO next time we're in this if
         * statement
         */
        sh4->last_inst_type = group;
    } else {
        /*
         * cash in on the dual-issue pipeline's "free" instruction and set
         * last_inst_type to SH4_GROUP_NONE so that the next instruction is
         * not free.
         */
        n_cycles = 0;
        sh4->last_inst_type = SH4_GROUP_NONE;
    }
    return n_cycles;
}

static inline void
sh4_do_exec_inst(Sh4 *sh4, inst_t inst, opcode_func_t op_fn) {
    Sh4OpArgs oa;
    oa.inst = inst;
    bool delayed_branch_tmp = sh4->delayed_branch;
    addr32_t delayed_branch_addr_tmp = sh4->delayed_branch_addr;

#ifdef DEEP_SYSCALL_TRACE
    deep_syscall_notify_jump(sh4->reg[SH4_REG_PC]);
#endif
    op_fn(sh4, oa);

#ifdef ENABLE_DEBUGGER
    if (!sh4->aborted_operation) {
        if (delayed_branch_tmp) {
            sh4->reg[SH4_REG_PC] = delayed_branch_addr_tmp;
            sh4->delayed_branch = false;

#ifdef DEEP_SYSCALL_TRACE
            deep_syscall_notify_jump(sh4->reg[SH4_REG_PC]);
#endif
        }
    } else {
        sh4->aborted_operation = false;
    }
#else
    if (delayed_branch_tmp) {
        sh4->reg[SH4_REG_PC] = delayed_branch_addr_tmp;
        sh4->delayed_branch = false;
    }
#endif
}

#endif
