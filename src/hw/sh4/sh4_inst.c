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

#include <assert.h>
#include <fenv.h>
#include <limits.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <stdbool.h>

#ifdef ENABLE_SH4_MMU
#include "sh4_mmu.h"
#endif

#include "error.h"
#include "dreamcast.h"
#include "sh4_ocache.h"
#include "sh4.h"
#include "sh4_tbl.h"
#include "sh4_excp.h"
#include "sh4_inst_def.h"
#include "log.h"

#ifdef ENABLE_DEBUGGER
#include "debugger.h"
#endif

#ifdef DEEP_SYSCALL_TRACE
#include "deep_syscall_trace.h"
#endif

#include "sh4_inst.h"

DEF_ERROR_STRING_ATTR(opcode_format)
DEF_ERROR_STRING_ATTR(opcode_name)
DEF_ERROR_INT_ATTR(instruction)
DEF_ERROR_INT_ATTR(instruction_mask)
DEF_ERROR_INT_ATTR(instruction_expect)
DEF_ERROR_U32_ATTR(fpscr)
DEF_ERROR_U32_ATTR(fpscr_expect)
DEF_ERROR_U32_ATTR(fpscr_mask)
DEF_ERROR_INT_ATTR(inst_bin)

#ifdef SH4_FPU_PEDANTIC
/*
 * set the FPU's invalid operation flag in FPSCR and maybe raise an exception
 *
 * dst is the destination register index for the operation.
 * It will be set to qNaN if the exceptions are disabled.
 */
static void sh4_fr_invalid(Sh4 *sh4, unsigned dst_reg);
#endif

static InstOpcode const* sh4_decode_inst_slow(inst_t inst);

static struct InstOpcode opcode_list[] = {
    // DONE
    // RTS
    { "0000000000001011", &sh4_inst_rts, true, SH4_GROUP_CO, 2 },

    // DONE
    // CLRMAC
    { "0000000000101000", &sh4_inst_clrmac, false, SH4_GROUP_CO, 1 },

    // DONE
    // CLRS
    { "0000000001001000", &sh4_inst_clrs, false, SH4_GROUP_CO, 1 },

    // DONE
    // CLRT
    { "0000000000001000", &sh4_inst_clrt, false, SH4_GROUP_MT, 1 },

    // DONE
    // LDTLB
    { "0000000000111000", &sh4_inst_ldtlb, false, SH4_GROUP_CO, 1 },

    // DONE
    // NOP
    { "0000000000001001", &sh4_inst_nop, false, SH4_GROUP_MT, 1 },

    // DONE
    // RTE
    { "0000000000101011", &sh4_inst_rte, false, SH4_GROUP_CO, 5 },

    // DONE
    // SETS
    { "0000000001011000", &sh4_inst_sets, false, SH4_GROUP_CO, 1 },

    // DONE
    // SETT
    { "0000000000011000", &sh4_inst_sett, false, SH4_GROUP_MT, 1 },

    // DONE
    // SLEEP
    { "0000000000011011", &sh4_inst_sleep, false, SH4_GROUP_CO, 4 },

    // DONE
    // FRCHG
    { "1111101111111101", &sh4_inst_frchg, false, SH4_GROUP_FE, 1 },

    // DONE
    // FSCHG
    { "1111001111111101", &sh4_inst_fschg, false, SH4_GROUP_FE, 1 },

    // DONE
    // MOVT Rn
    { "0000nnnn00101001", &sh4_inst_unary_movt_gen, false,
      SH4_GROUP_EX, 1 },

    // DONE
    // CMP/PZ
    { "0100nnnn00010001", &sh4_inst_unary_cmppz_gen, false,
      SH4_GROUP_MT, 1 },

    // DONE
    // CMP/PL
    { "0100nnnn00010101", &sh4_inst_unary_cmppl_gen, false,
      SH4_GROUP_MT, 1 },

    // DONE
    // DT
    { "0100nnnn00010000", &sh4_inst_unary_dt_gen, false,
      SH4_GROUP_EX, 1 },

    // DONE
    // ROTL Rn
    { "0100nnnn00000100", &sh4_inst_unary_rotl_gen, false,
      SH4_GROUP_EX, 1 },

    // DONE
    // ROTR Rn
    { "0100nnnn00000101", &sh4_inst_unary_rotr_gen, false,
      SH4_GROUP_EX, 1 },

    // DONE
    // ROTCL Rn
    { "0100nnnn00100100", &sh4_inst_unary_rotcl_gen, false,
      SH4_GROUP_EX, 1 },

    // DONE
    // ROTCR Rn
    { "0100nnnn00100101", &sh4_inst_unary_rotcr_gen, false,
      SH4_GROUP_EX, 1 },

    // DONE
    // SHAL Rn
    { "0100nnnn00100000", &sh4_inst_unary_shal_gen, false,
      SH4_GROUP_EX, 1 },

    // DONE
    // SHAR Rn
    { "0100nnnn00100001", &sh4_inst_unary_shar_gen, false,
      SH4_GROUP_EX, 1 },

    // DONE
    // SHLL Rn
    { "0100nnnn00000000", &sh4_inst_unary_shll_gen, false,
      SH4_GROUP_EX, 1 },

    // DONE
    // SHLR Rn
    { "0100nnnn00000001", &sh4_inst_unary_shlr_gen, false,
      SH4_GROUP_EX, 1 },

    // DONE
    // SHLL2 Rn
    { "0100nnnn00001000", &sh4_inst_unary_shll2_gen, false,
      SH4_GROUP_EX, 1 },

    // DONE
    // SHLR2 Rn
    { "0100nnnn00001001", &sh4_inst_unary_shlr2_gen, false,
      SH4_GROUP_EX, 1 },

    // DONE
    // SHLL8 Rn
    { "0100nnnn00011000", &sh4_inst_unary_shll8_gen, false,
      SH4_GROUP_EX, 1 },

    // DONE
    // SHLR8 Rn
    { "0100nnnn00011001", &sh4_inst_unary_shlr8_gen, false,
      SH4_GROUP_EX, 1 },

    // DONE
    // SHLL16 Rn
    { "0100nnnn00101000", &sh4_inst_unary_shll16_gen, false,
      SH4_GROUP_EX, 1 },

    // DONE
    // SHLR16 Rn
    { "0100nnnn00101001", &sh4_inst_unary_shlr16_gen, false,
      SH4_GROUP_EX, 1 },

    // DONE
    // BRAF Rn
    { "0000nnnn00100011", &sh4_inst_unary_braf_gen, true,
      SH4_GROUP_CO, 2 },

    // DONE
    // BSRF Rn
    { "0000nnnn00000011", &sh4_inst_unary_bsrf_gen, true,
      SH4_GROUP_CO, 2 },

    // CMP/EQ #imm, R0
    { "10001000iiiiiiii", &sh4_inst_binary_cmpeq_imm_r0, false,
      SH4_GROUP_MT, 1 },

    // AND.B #imm, @(R0, GBR)
    { "11001101iiiiiiii", &sh4_inst_binary_andb_imm_r0_gbr, false,
      SH4_GROUP_CO, 4 },

    // AND #imm, R0
    { "11001001iiiiiiii", &sh4_inst_binary_and_imm_r0, false,
      SH4_GROUP_EX, 1 },

    // OR.B #imm, @(R0, GBR)
    { "11001111iiiiiiii", &sh4_inst_binary_orb_imm_r0_gbr, false,
      SH4_GROUP_CO, 4 },

    // OR #imm, R0
    { "11001011iiiiiiii", &sh4_inst_binary_or_imm_r0, false,
      SH4_GROUP_EX, 1 },

    // TST #imm, R0
    { "11001000iiiiiiii", &sh4_inst_binary_tst_imm_r0, false,
      SH4_GROUP_MT, 1 },

    // TST.B #imm, @(R0, GBR)
    { "11001100iiiiiiii", &sh4_inst_binary_tstb_imm_r0_gbr, false,
      SH4_GROUP_CO, 3 },

    // XOR #imm, R0
    { "11001010iiiiiiii", &sh4_inst_binary_xor_imm_r0, false,
      SH4_GROUP_EX, 1 },

    // XOR.B #imm, @(R0, GBR)
    { "11001110iiiiiiii", &sh4_inst_binary_xorb_imm_r0_gbr, false,
      SH4_GROUP_CO, 4 },

    // BF label
    { "10001011dddddddd", &sh4_inst_unary_bf_disp, true,
      SH4_GROUP_BR, 1 },

    // BF/S label
    { "10001111dddddddd", &sh4_inst_unary_bfs_disp, true,
      SH4_GROUP_BR, 1 },

    // BT label
    { "10001001dddddddd", &sh4_inst_unary_bt_disp, true,
      SH4_GROUP_BR, 1 },

    // BT/S label
    { "10001101dddddddd", &sh4_inst_unary_bts_disp, true,
      SH4_GROUP_BR, 1 },

    // BRA label
    { "1010dddddddddddd", &sh4_inst_unary_bra_disp, true,
      SH4_GROUP_BR, 1 },

    // BSR label
    { "1011dddddddddddd", &sh4_inst_unary_bsr_disp, true,
      SH4_GROUP_BR, 1 },

    // TRAPA #immed
    { "11000011iiiiiiii", &sh4_inst_unary_trapa_disp, false,
      SH4_GROUP_CO, 7 },

    // DONE
    // TAS.B @Rn
    { "0100nnnn00011011", &sh4_inst_unary_tasb_gen, false,
      SH4_GROUP_CO, 5 },

    // DONE
    // OCBI @Rn
    { "0000nnnn10010011", &sh4_inst_unary_ocbi_indgen, false,
      SH4_GROUP_LS, 1 },

    // DONE
    // OCBP @Rn
    { "0000nnnn10100011", &sh4_inst_unary_ocbp_indgen, false,
      SH4_GROUP_LS, 1 },

    // DONE
    // OCBWB @Rn
    { "0000nnnn10110011", &sh4_inst_unary_ocbwb_indgen, false,
      SH4_GROUP_LS, 1 },

    // DONE
    // PREF @Rn
    { "0000nnnn10000011", &sh4_inst_unary_pref_indgen, false,
      SH4_GROUP_LS, 1 },

    // DONE
    // JMP @Rn
    { "0100nnnn00101011", &sh4_inst_unary_jmp_indgen, true,
      SH4_GROUP_CO, 2 },

    // DONE
    // JSR @Rn
    { "0100nnnn00001011", &sh4_inst_unary_jsr_indgen, true,
      SH4_GROUP_CO, 2 },

    // DONE
    // LDC Rm, SR
    { "0100mmmm00001110", &sh4_inst_binary_ldc_gen_sr, false,
      SH4_GROUP_CO, 4 },

    // DONE
    // LDC Rm, GBR
    { "0100mmmm00011110", &sh4_inst_binary_ldc_gen_gbr, false,
      SH4_GROUP_CO, 3 },

    // DONE
    // LDC Rm, VBR
    { "0100mmmm00101110", &sh4_inst_binary_ldc_gen_vbr, false,
      SH4_GROUP_CO, 1 },

    // DONE
    // LDC Rm, SSR
    { "0100mmmm00111110", &sh4_inst_binary_ldc_gen_ssr, false,
      SH4_GROUP_CO, 1 },

    // DONE
    // LDC Rm, SPC
    { "0100mmmm01001110", &sh4_inst_binary_ldc_gen_spc, false,
      SH4_GROUP_CO, 1 },

    // DONE
    // LDC Rm, DBR
    { "0100mmmm11111010", &sh4_inst_binary_ldc_gen_dbr, false,
      SH4_GROUP_CO, 1 },

    // DONE
    // STC SR, Rn
    { "0000nnnn00000010", &sh4_inst_binary_stc_sr_gen, false,
      SH4_GROUP_CO, 2 },

    // DONE
    // STC GBR, Rn
    { "0000nnnn00010010", &sh4_inst_binary_stc_gbr_gen, false,
      SH4_GROUP_CO, 2 },

    // DONE
    // STC VBR, Rn
    { "0000nnnn00100010", &sh4_inst_binary_stc_vbr_gen, false,
      SH4_GROUP_CO, 2 },

    // DONE
    // STC SSR, Rn
    { "0000nnnn00110010", &sh4_inst_binary_stc_ssr_gen, false,
      SH4_GROUP_CO, 2 },

    // DONE
    // STC SPC, Rn
    { "0000nnnn01000010", &sh4_inst_binary_stc_spc_gen, false,
      SH4_GROUP_CO, 2 },

    // DONE
    // STC SGR, Rn
    { "0000nnnn00111010", &sh4_inst_binary_stc_sgr_gen, false,
      SH4_GROUP_CO, 3 },

    // DONE
    // STC DBR, Rn
    { "0000nnnn11111010", &sh4_inst_binary_stc_dbr_gen, false,
      SH4_GROUP_CO, 2 },

    // DONE
    // LDC.L @Rm+, SR
    { "0100mmmm00000111", &sh4_inst_binary_ldcl_indgeninc_sr, false,
      SH4_GROUP_CO, 4 },

    // DONE
    // LDC.L @Rm+, GBR
    { "0100mmmm00010111", &sh4_inst_binary_ldcl_indgeninc_gbr, false,
      SH4_GROUP_CO, 3 },

    // DONE
    // LDC.L @Rm+, VBR
    { "0100mmmm00100111", &sh4_inst_binary_ldcl_indgeninc_vbr, false,
      SH4_GROUP_CO, 1 },

    // DONE
    // LDC.L @Rm+, SSR
    { "0100mmmm00110111", &sh4_inst_binary_ldcl_indgenic_ssr, false,
      SH4_GROUP_CO, 1 },

    // DONE
    // LDC.L @Rm+, SPC
    { "0100mmmm01000111", &sh4_inst_binary_ldcl_indgeninc_spc, false,
      SH4_GROUP_CO, 1 },

    // DONE
    // LDC.L @Rm+, DBR
    { "0100mmmm11110110", &sh4_inst_binary_ldcl_indgeninc_dbr, false,
      SH4_GROUP_CO, 1 },

    // DONE
    // STC.L SR, @-Rn
    { "0100nnnn00000011", &sh4_inst_binary_stcl_sr_inddecgen, false,
      SH4_GROUP_CO, 2 },

    // DONE
    // STC.L GBR, @-Rn
    { "0100nnnn00010011", &sh4_inst_binary_stcl_gbr_inddecgen, false,
      SH4_GROUP_CO, 2 },

    // DONE
    // STC.L VBR, @-Rn
    { "0100nnnn00100011", &sh4_inst_binary_stcl_vbr_inddecgen, false,
      SH4_GROUP_CO, 2 },

    // DONE
    // STC.L SSR, @-Rn
    { "0100nnnn00110011", &sh4_inst_binary_stcl_ssr_inddecgen, false,
      SH4_GROUP_CO, 2 },

    // DONE
    // STC.L SPC, @-Rn
    { "0100nnnn01000011", &sh4_inst_binary_stcl_spc_inddecgen, false,
      SH4_GROUP_CO, 2 },

    // DONE
    // STC.L SGR, @-Rn
    { "0100nnnn00110010", &sh4_inst_binary_stcl_sgr_inddecgen, false,
      SH4_GROUP_CO, 3 },

    // DONE
    // STC.L DBR, @-Rn
    { "0100nnnn11110010", &sh4_inst_binary_stcl_dbr_inddecgen, false,
      SH4_GROUP_CO, 2 },

    // MOV #imm, Rn
    { "1110nnnniiiiiiii", &sh4_inst_binary_mov_imm_gen, false,
      SH4_GROUP_EX, 1 },

    // ADD #imm, Rn
    { "0111nnnniiiiiiii", &sh4_inst_binary_add_imm_gen, false,
      SH4_GROUP_EX, 1 },

    // MOV.W @(disp, PC), Rn
    { "1001nnnndddddddd", &sh4_inst_binary_movw_binind_disp_pc_gen,
      false, SH4_GROUP_LS, 1 },

    // MOV.L @(disp, PC), Rn
    { "1101nnnndddddddd", &sh4_inst_binary_movl_binind_disp_pc_gen,
      false, SH4_GROUP_LS, 1 },

    // DONE
    // MOV Rm, Rn
    { "0110nnnnmmmm0011", &sh4_inst_binary_mov_gen_gen, false,
      SH4_GROUP_MT, 1 },

    // DONE
    // SWAP.B Rm, Rn
    { "0110nnnnmmmm1000", &sh4_inst_binary_swapb_gen_gen, false,
      SH4_GROUP_EX, 1 },

    // DONE
    // SWAP.W Rm, Rn
    { "0110nnnnmmmm1001", &sh4_inst_binary_swapw_gen_gen, false,
      SH4_GROUP_EX, 1 },

    // DONE
    // XTRCT Rm, Rn
    { "0010nnnnmmmm1101", &sh4_inst_binary_xtrct_gen_gen, false,
      SH4_GROUP_EX, 1 },

    // DONE
    // ADD Rm, Rn
    { "0011nnnnmmmm1100", &sh4_inst_binary_add_gen_gen, false,
      SH4_GROUP_EX, 1 },

    // DONE
    // ADDC Rm, Rn
    { "0011nnnnmmmm1110", &sh4_inst_binary_addc_gen_gen, false,
      SH4_GROUP_EX, 1 },

    // DONE
    // ADDV Rm, Rn
    { "0011nnnnmmmm1111", &sh4_inst_binary_addv_gen_gen, false,
      SH4_GROUP_EX, 1 },

    // DONE
    // CMP/EQ Rm, Rn
    { "0011nnnnmmmm0000", &sh4_inst_binary_cmpeq_gen_gen, false,
      SH4_GROUP_MT, 1 },

    // DONE
    // CMP/HS Rm, Rn
    { "0011nnnnmmmm0010", &sh4_inst_binary_cmphs_gen_gen, false,
      SH4_GROUP_MT, 1 },

    // DONE
    // CMP/GE Rm, Rn
    { "0011nnnnmmmm0011", &sh4_inst_binary_cmpge_gen_gen, false,
      SH4_GROUP_MT, 1 },

    // DONE
    // CMP/HI Rm, Rn
    { "0011nnnnmmmm0110", &sh4_inst_binary_cmphi_gen_gen, false,
      SH4_GROUP_MT, 1 },

    // DONE
    // CMP/GT Rm, Rn
    { "0011nnnnmmmm0111", &sh4_inst_binary_cmpgt_gen_gen, false,
      SH4_GROUP_MT, 1 },

    // DONE
    // CMP/STR Rm, Rn
    { "0010nnnnmmmm1100", &sh4_inst_binary_cmpstr_gen_gen, false,
      SH4_GROUP_MT, 1 },

    // DONE
    // DIV1 Rm, Rn
    { "0011nnnnmmmm0100", &sh4_inst_binary_div1_gen_gen, false,
      SH4_GROUP_EX, 1 },

    // DONE
    // DIV0S Rm, Rn
    { "0010nnnnmmmm0111", &sh4_inst_binary_div0s_gen_gen, false,
      SH4_GROUP_EX, 1 },

    // DIV0U
    { "0000000000011001", &sh4_inst_noarg_div0u, false, SH4_GROUP_EX, 1 },

    // DONE
    // DMULS.L Rm, Rn
    { "0011nnnnmmmm1101", &sh4_inst_binary_dmulsl_gen_gen, false,
      SH4_GROUP_CO, 2 },

    // DONE
    // DMULU.L Rm, Rn
    { "0011nnnnmmmm0101", &sh4_inst_binary_dmulul_gen_gen, false,
      SH4_GROUP_CO, 2 },

    // DONE
    // EXTS.B Rm, Rn
    { "0110nnnnmmmm1110", &sh4_inst_binary_extsb_gen_gen, false,
      SH4_GROUP_EX, 1 },

    // DONE
    // EXTS.W Rm, Rn
    { "0110nnnnmmmm1111", &sh4_inst_binary_extsw_gen_gen, false,
      SH4_GROUP_EX, 1 },

    // DONE
    // EXTU.B Rm, Rn
    { "0110nnnnmmmm1100", &sh4_inst_binary_extub_gen_gen, false,
      SH4_GROUP_EX, 1 },

    // DONE
    // EXTU.W Rm, Rn
    { "0110nnnnmmmm1101", &sh4_inst_binary_extuw_gen_gen, false,
      SH4_GROUP_EX, 1 },

    // DONE
    // MUL.L Rm, Rn
    { "0000nnnnmmmm0111", &sh4_inst_binary_mull_gen_gen, false,
      SH4_GROUP_CO, 2 },

    // DONE
    // MULS.W Rm, Rn
    { "0010nnnnmmmm1111", &sh4_inst_binary_mulsw_gen_gen, false,
      SH4_GROUP_CO, 2 },

    // DONE
    // MULU.W Rm, Rn
    { "0010nnnnmmmm1110", &sh4_inst_binary_muluw_gen_gen, false,
      SH4_GROUP_CO, 2 },

    // DONE
    // NEG Rm, Rn
    { "0110nnnnmmmm1011", &sh4_inst_binary_neg_gen_gen, false,
      SH4_GROUP_EX, 1 },

    // DONE
    // NEGC Rm, Rn
    { "0110nnnnmmmm1010", &sh4_inst_binary_negc_gen_gen, false,
      SH4_GROUP_EX, 1 },

    // DONE
    // SUB Rm, Rn
    { "0011nnnnmmmm1000", &sh4_inst_binary_sub_gen_gen, false,
      SH4_GROUP_EX, 1 },

    // DONE
    // SUBC Rm, Rn
    { "0011nnnnmmmm1010", &sh4_inst_binary_subc_gen_gen, false,
      SH4_GROUP_EX, 1 },

    // DONE
    // SUBV Rm, Rn
    { "0011nnnnmmmm1011", &sh4_inst_binary_subv_gen_gen, false,
      SH4_GROUP_EX, 1 },

    // DONE
    // AND Rm, Rn
    { "0010nnnnmmmm1001", &sh4_inst_binary_and_gen_gen, false,
      SH4_GROUP_EX, 1 },

    // DONE
    // NOT Rm, Rn
    { "0110nnnnmmmm0111", &sh4_inst_binary_not_gen_gen, false,
      SH4_GROUP_EX, 1 },

    // DONE
    // OR Rm, Rn
    { "0010nnnnmmmm1011", &sh4_inst_binary_or_gen_gen, false,
      SH4_GROUP_EX, 1 },

    // DONE
    // TST Rm, Rn
    { "0010nnnnmmmm1000", &sh4_inst_binary_tst_gen_gen, false,
      SH4_GROUP_MT, 1 },

    // DONE
    // XOR Rm, Rn
    { "0010nnnnmmmm1010", &sh4_inst_binary_xor_gen_gen, false,
      SH4_GROUP_EX, 1 },

    // DONE
    // SHAD Rm, Rn
    { "0100nnnnmmmm1100", &sh4_inst_binary_shad_gen_gen, false,
      SH4_GROUP_EX, 1 },

    // DONE
    // SHLD Rm, Rn
    { "0100nnnnmmmm1101", &sh4_inst_binary_shld_gen_gen, false,
      SH4_GROUP_EX, 1 },

    // LDC Rm, Rn_BANK
    { "0100mmmm1nnn1110", &sh4_inst_binary_ldc_gen_bank, false,
      SH4_GROUP_CO, 1 },

    // LDC.L @Rm+, Rn_BANK
    { "0100mmmm1nnn0111", &sh4_inst_binary_ldcl_indgeninc_bank, false,
      SH4_GROUP_CO, 1 },

    // STC Rm_BANK, Rn
    { "0000nnnn1mmm0010", &sh4_inst_binary_stc_bank_gen, false,
      SH4_GROUP_CO, 2 },

    // STC.L Rm_BANK, @-Rn
    { "0100nnnn1mmm0011", &sh4_inst_binary_stcl_bank_inddecgen, false,
      SH4_GROUP_CO, 2 },

    // DONE
    // LDS Rm, MACH
    { "0100mmmm00001010", &sh4_inst_binary_lds_gen_mach, false,
      SH4_GROUP_CO, 1 },

    // DONE
    // LDS Rm, MACL
    { "0100mmmm00011010", &sh4_inst_binary_lds_gen_macl, false,
      SH4_GROUP_CO, 1 },

    // DONE
    // STS MACH, Rn
    { "0000nnnn00001010", &sh4_inst_binary_sts_mach_gen, false,
      SH4_GROUP_CO, 1 },

    // DONE
    // STS MACL, Rn
    { "0000nnnn00011010", &sh4_inst_binary_sts_macl_gen, false,
      SH4_GROUP_CO, 1 },

    // DONE
    // LDS Rm, PR
    { "0100mmmm00101010", &sh4_inst_binary_lds_gen_pr, false,
      SH4_GROUP_CO, 2 },

    // DONE
    // STS PR, Rn
    { "0000nnnn00101010", &sh4_inst_binary_sts_pr_gen, false,
      SH4_GROUP_CO, 2 },

    // DONE
    // LDS.L @Rm+, MACH
    { "0100mmmm00000110", &sh4_inst_binary_ldsl_indgeninc_mach, false,
      SH4_GROUP_CO, 1 },

    // DONE
    // LDS.L @Rm+, MACL
    { "0100mmmm00010110", &sh4_inst_binary_ldsl_indgeninc_macl, false,
      SH4_GROUP_CO, 1 },

    // DONE
    // STS.L MACH, @-Rn
    { "0100mmmm00000010", &sh4_inst_binary_stsl_mach_inddecgen, false,
      SH4_GROUP_CO, 1 },

    // DONE
    // STS.L MACL, @-Rn
    { "0100mmmm00010010", &sh4_inst_binary_stsl_macl_inddecgen, false,
      SH4_GROUP_CO, 1 },

    // DONE
    // LDS.L @Rm+, PR
    { "0100mmmm00100110", &sh4_inst_binary_ldsl_indgeninc_pr, false,
      SH4_GROUP_CO, 2 },

    // DONE
    // STS.L PR, @-Rn
    { "0100nnnn00100010", &sh4_inst_binary_stsl_pr_inddecgen, false,
      SH4_GROUP_CO, 2 },

    // DONE
    // MOV.B Rm, @Rn
    { "0010nnnnmmmm0000", &sh4_inst_binary_movb_gen_indgen, false,
      SH4_GROUP_LS, 1 },

    // DONE
    // MOV.W Rm, @Rn
    { "0010nnnnmmmm0001", &sh4_inst_binary_movw_gen_indgen, false,
      SH4_GROUP_LS, 1 },

    // DONE
    // MOV.L Rm, @Rn
    { "0010nnnnmmmm0010", &sh4_inst_binary_movl_gen_indgen, false,
      SH4_GROUP_LS, 1 },

    // DONE
    // MOV.B @Rm, Rn
    { "0110nnnnmmmm0000", &sh4_inst_binary_movb_indgen_gen, false,
      SH4_GROUP_LS, 1 },

    // DONE
    // MOV.W @Rm, Rn
    { "0110nnnnmmmm0001", &sh4_inst_binary_movw_indgen_gen, false,
      SH4_GROUP_LS, 1 },

    // DONE
    // MOV.L @Rm, Rn
    { "0110nnnnmmmm0010", &sh4_inst_binary_movl_indgen_gen, false,
      SH4_GROUP_LS, 1 },

    // DONE
    // MOV.B Rm, @-Rn
    { "0010nnnnmmmm0100", &sh4_inst_binary_movb_gen_inddecgen, false,
      SH4_GROUP_LS, 1 },

    // DONE
    // MOV.W Rm, @-Rn
    { "0010nnnnmmmm0101", &sh4_inst_binary_movw_gen_inddecgen, false,
      SH4_GROUP_LS, 1 },

    // DONE
    // MOV.L Rm, @-Rn
    { "0010nnnnmmmm0110", &sh4_inst_binary_movl_gen_inddecgen, false,
      SH4_GROUP_LS, 1 },

    // DONE
    // MOV.B @Rm+, Rn
    { "0110nnnnmmmm0100", &sh4_inst_binary_movb_indgeninc_gen, false,
      SH4_GROUP_LS, 1 },

    // DONE
    // MOV.W @Rm+, Rn
    { "0110nnnnmmmm0101", &sh4_inst_binary_movw_indgeninc_gen, false,
      SH4_GROUP_LS, 1 },

    // DONE
    // MOV.L @Rm+, Rn
    { "0110nnnnmmmm0110", &sh4_inst_binary_movl_indgeninc_gen, false,
      SH4_GROUP_LS, 1 },

    // DONE
    // MAC.L @Rm+, @Rn+
    { "0000nnnnmmmm1111", &sh4_inst_binary_macl_indgeninc_indgeninc,
      false, SH4_GROUP_CO, 2 },

    // DONE
    // MAC.W @Rm+, @Rn+
    { "0100nnnnmmmm1111", &sh4_inst_binary_macw_indgeninc_indgeninc,
      false, SH4_GROUP_CO, 2 },

    // MOV.B R0, @(disp, Rn)
    { "10000000nnnndddd", &sh4_inst_binary_movb_r0_binind_disp_gen,
      false, SH4_GROUP_LS, 1 },

    // MOV.W R0, @(disp, Rn)
    { "10000001nnnndddd", &sh4_inst_binary_movw_r0_binind_disp_gen,
      false, SH4_GROUP_LS, 1 },

    // MOV.L Rm, @(disp, Rn)
    { "0001nnnnmmmmdddd", &sh4_inst_binary_movl_gen_binind_disp_gen,
      false, SH4_GROUP_LS, 1 },

    // MOV.B @(disp, Rm), R0
    { "10000100mmmmdddd", &sh4_inst_binary_movb_binind_disp_gen_r0,
      false, SH4_GROUP_LS, 1 },

    // MOV.W @(disp, Rm), R0
    { "10000101mmmmdddd", &sh4_inst_binary_movw_binind_disp_gen_r0,
      false, SH4_GROUP_LS, 1 },

    // MOV.L @(disp, Rm), Rn
    { "0101nnnnmmmmdddd", &sh4_inst_binary_movl_binind_disp_gen_gen,
      false, SH4_GROUP_LS, 1 },

    // DONE
    // MOV.B Rm, @(R0, Rn)
    { "0000nnnnmmmm0100", &sh4_inst_binary_movb_gen_binind_r0_gen,
      false, SH4_GROUP_LS, 1 },

    // DONE
    // MOV.W Rm, @(R0, Rn)
    { "0000nnnnmmmm0101", &sh4_inst_binary_movw_gen_binind_r0_gen,
      false, SH4_GROUP_LS, 1 },

    // DONE
    // MOV.L Rm, @(R0, Rn)
    { "0000nnnnmmmm0110", &sh4_inst_binary_movl_gen_binind_r0_gen,
      false, SH4_GROUP_LS, 1 },

    // DONE
    // MOV.B @(R0, Rm), Rn
    { "0000nnnnmmmm1100", &sh4_inst_binary_movb_binind_r0_gen_gen,
      false, SH4_GROUP_LS, 1 },

    // DONE
    // MOV.W @(R0, Rm), Rn
    { "0000nnnnmmmm1101", &sh4_inst_binary_movw_binind_r0_gen_gen,
      false, SH4_GROUP_LS, 1 },

    // DONE
    // MOV.L @(R0, Rm), Rn
    { "0000nnnnmmmm1110", &sh4_inst_binary_movl_binind_r0_gen_gen,
      false, SH4_GROUP_LS, 1 },

    // MOV.B R0, @(disp, GBR)
    { "11000000dddddddd", &sh4_inst_binary_movb_r0_binind_disp_gbr,
      false, SH4_GROUP_LS, 1 },

    // MOV.W R0, @(disp, GBR)
    { "11000001dddddddd", &sh4_inst_binary_movw_r0_binind_disp_gbr,
      false, SH4_GROUP_LS, 1 },

    // MOV.L R0, @(disp, GBR)
    { "11000010dddddddd", &sh4_inst_binary_movl_r0_binind_disp_gbr,
      false, SH4_GROUP_LS, 1 },

    // MOV.B @(disp, GBR), R0
    { "11000100dddddddd", &sh4_inst_binary_movb_binind_disp_gbr_r0,
      false, SH4_GROUP_LS, 1 },

    // MOV.W @(disp, GBR), R0
    { "11000101dddddddd", &sh4_inst_binary_movw_binind_disp_gbr_r0,
      false, SH4_GROUP_LS, 1 },

    // MOV.L @(disp, GBR), R0
    { "11000110dddddddd", &sh4_inst_binary_movl_binind_disp_gbr_r0,
      false, SH4_GROUP_LS, 1 },

    // MOVA @(disp, PC), R0
    { "11000111dddddddd", &sh4_inst_binary_mova_binind_disp_pc_r0,
      false, SH4_GROUP_EX, 1 },

    // MOVCA.L R0, @Rn
    { "0000nnnn11000011", &sh4_inst_binary_movcal_r0_indgen,
      false, SH4_GROUP_LS, 1 },

    // DONE
    // FLDI0 FRn
    { "1111nnnn10001101", FPU_HANDLER(fldi0),
      false, SH4_GROUP_LS, 1 },

    // DONE
    // FLDI1 Frn
    { "1111nnnn10011101", FPU_HANDLER(fldi1),
      false, SH4_GROUP_LS, 1 },

    // DONE
    // FMOV FRm, FRn
    // 1111nnnnmmmm1100
    // FMOV DRm, DRn
    // 1111nnn0mmm01100
    // FMOV XDm, DRn
    // 1111nnn0mmm11100
    // FMOV DRm, XDn
    // 1111nnn1mmm01100
    // FMOV XDm, XDn
    // 1111nnn1mmm11100
    { "1111nnnnmmmm1100", FPU_HANDLER(fmov_gen),
      false, SH4_GROUP_LS, 1 },

    // DONE
    // FMOV.S @Rm, FRn
    // 1111nnnnmmmm1000
    // FMOV @Rm, DRn
    // 1111nnn0mmmm1000
    // FMOV @Rm, XDn
    // 1111nnn1mmmm1000
    { "1111nnnnmmmm1000", FPU_HANDLER(fmovs_ind_gen),
      false, SH4_GROUP_LS, 1 },

    // DONE
    // FMOV.S @(R0, Rm), FRn
    // 1111nnnnmmmm0110
    // FMOV @(R0, Rm), DRn
    // 1111nnn0mmmm0110
    // FMOV @(R0, Rm), XDn
    // 1111nnn1mmmm0110
    { "1111nnnnmmmm0110", FPU_HANDLER(fmov_binind_r0_gen_fpu),
      false, SH4_GROUP_LS, 1 },

    // DONE
    // FMOV.S @Rm+, FRn
    // 1111nnnnmmmm1001
    // FMOV @Rm+, DRn
    // 1111nnn0mmmm1001
    // FMOV @Rm+, XDn
    // 1111nnn1mmmm1001
    { "1111nnnnmmmm1001", FPU_HANDLER(fmov_indgeninc_fpu),
      false, SH4_GROUP_LS, 1 },

    // DONE
    // FMOV.S FRm, @Rn
    // 1111nnnnmmmm1010
    // FMOV DRm, @Rn
    // 1111nnnnmmm01010
    // FMOV XDm, @Rn
    // 1111nnnnmmm11010
    { "1111nnnnmmmm1010", FPU_HANDLER(fmov_fpu_indgen),
      false, SH4_GROUP_LS, 1 },

    // DONE
    // FMOV.S FRm, @-Rn
    // 1111nnnnmmmm1011
    // FMOV DRm, @-Rn
    // 1111nnnnmmm01011
    // FMOV XDm, @-Rn
    // 1111nnnnmmm11011
    { "1111nnnnmmmm1011", FPU_HANDLER(fmov_fpu_inddecgen),
      false, SH4_GROUP_LS, 1 },

    // DONE
    // FMOV.S FRm, @(R0, Rn)
    // 1111nnnnmmmm0111
    // FMOV DRm, @(R0, Rn)
    // 1111nnnnmmm00111
    // FMOV XDm, @(R0, Rn)
    // 1111nnnnmmm10111
    { "1111nnnnmmmm0111", FPU_HANDLER(fmov_fpu_binind_r0_gen),
      false, SH4_GROUP_LS, 1 },

    // FLDS FRm, FPUL
    // XXX Should this check the SZ or PR bits of FPSCR ?
    { "1111mmmm00011101", &sh4_inst_binary_flds_fr_fpul, false,
      SH4_GROUP_LS, 1 },

    // FSTS FPUL, FRn
    // XXX Should this check the SZ or PR bits of FPSCR ?
    { "1111nnnn00001101", &sh4_inst_binary_fsts_fpul_fr, false,
      SH4_GROUP_LS, 1 },

    // FABS FRn
    // 1111nnnn01011101
    // FABS DRn
    // 1111nnn001011101
    { "1111nnnn01011101", FPU_HANDLER(fabs_fpu), false,
      SH4_GROUP_LS, 1 },

    // DONE
    // FADD FRm, FRn
    // 1111nnnnmmmm0000
    // FADD DRm, DRn
    // 1111nnn0mmm00000
    { "1111nnnnmmmm0000", FPU_HANDLER(fadd_fpu),
      false, SH4_GROUP_FE, 1 },

    // DONE
    // FCMP/EQ FRm, FRn
    // 1111nnnnmmmm0100
    // FCMP/EQ DRm, DRn
    // 1111nnn0mmm00100
    { "1111nnnnmmmm0100", FPU_HANDLER(fcmpeq_fpu), false, SH4_GROUP_FE, 1 },

    // DONE
    // FCMP/GT FRm, FRn
    // 1111nnnnmmmm0101
    // FCMP/GT DRm, DRn
    // 1111nnn0mmm00101
    { "1111nnnnmmmm0101", FPU_HANDLER(fcmpgt_fpu), false, SH4_GROUP_FE, 1 },

    // DONE
    // FDIV FRm, FRn
    // 1111nnnnmmmm0011
    // FDIV DRm, DRn
    // 1111nnn0mmm00011
    { "1111nnnnmmmm0011", FPU_HANDLER(fdiv_fpu), false, SH4_GROUP_FE, 1 },

    // FLOAT FPUL, FRn
    // 1111nnnn00101101
    // FLOAT FPUL, DRn
    // 1111nnn000101101
    { "1111nnnn00101101", FPU_HANDLER(float_fpu), false, SH4_GROUP_FE, 1 },

    // DONE
    // FMAC FR0, FRm, FRn
    // 1111nnnnmmmm1110
    { "1111nnnnmmmm1110", FPU_HANDLER(fmac_fpu), false, SH4_GROUP_FE, 1 },

    // DONE
    // FMUL FRm, FRn
    // 1111nnnnmmmm0010
    // FMUL DRm, DRn
    // 1111nnn0mmm00010
    { "1111nnnnmmmm0010", FPU_HANDLER(fmul_fpu), false, SH4_GROUP_FE, 1 },

    // FNEG FRn
    // 1111nnnn01001101
    // FNEG DRn
    // 1111nnn001001101
    { "1111nnnn01001101", FPU_HANDLER(fneg_fpu), false, SH4_GROUP_LS, 1 },

    // FSQRT FRn
    // 1111nnnn01101101
    // FSQRT DRn
    // 1111nnn001101101
    { "1111nnnn01101101", FPU_HANDLER(fsqrt_fpu), false, SH4_GROUP_FE, 1 },

    // DONE
    // FSUB FRm, FRn
    // 1111nnnnmmmm0001
    // FSUB DRm, DRn
    // 1111nnn0mmm00001
    { "1111nnnnmmmm0001", FPU_HANDLER(fsub_fpu), false, SH4_GROUP_FE, 1 },

    // FTRC FRm, FPUL
    // 1111mmmm00111101
    // FTRC DRm, FPUL
    // 1111mmm000111101
    { "1111mmmm00111101", FPU_HANDLER(ftrc_fpu), false, SH4_GROUP_FE, 1 },

    // FCNVDS DRm, FPUL
    // 1111mmm010111101
    { "1111mmm010111101", FPU_HANDLER(fcnvds_fpu), false, SH4_GROUP_FE, 1 },

    // FCNVSD FPUL, DRn
    // 1111nnn010101101
    { "1111nnn010101101", FPU_HANDLER(fcnvsd_fpu), false, SH4_GROUP_FE, 1 },

    // LDS Rm, FPSCR
    { "0100mmmm01101010", &sh4_inst_binary_lds_gen_fpscr, false,
      SH4_GROUP_CO, 1 },

    // LDS Rm, FPUL
    { "0100mmmm01011010", &sh4_inst_binary_gen_fpul, false,
      SH4_GROUP_LS, 1 },

    // LDS.L @Rm+, FPSCR
    { "0100mmmm01100110", &sh4_inst_binary_ldsl_indgeninc_fpscr, false,
      SH4_GROUP_CO, 1 },

    // LDS.L @Rm+, FPUL
    { "0100mmmm01010110", &sh4_inst_binary_ldsl_indgeninc_fpul, false,
      SH4_GROUP_CO, 1 },

    // STS FPSCR, Rn
    { "0000nnnn01101010", &sh4_inst_binary_sts_fpscr_gen, false,
      SH4_GROUP_CO, 1 },

    // STS FPUL, Rn
    { "0000nnnn01011010", &sh4_inst_binary_sts_fpul_gen, false,
      SH4_GROUP_LS, 1 },

    // STS.L FPSCR, @-Rn
    { "0100nnnn01100010", &sh4_inst_binary_stsl_fpscr_inddecgen, false,
      SH4_GROUP_CO, 1 },

    // STS.L FPUL, @-Rn
    { "0100nnnn01010010", &sh4_inst_binary_stsl_fpul_inddecgen, false,
      SH4_GROUP_CO, 1 },

    // FIPR FVm, FVn - vector dot product
    { "1111nnmm11101101", &sh4_inst_binary_fipr_fv_fv, false,
      SH4_GROUP_FE, 1 },

    // FTRV XMTRX, FVn - multiple vector by matrix
    { "1111nn0111111101", &sh4_inst_binary_fitrv_mxtrx_fv, false,
      SH4_GROUP_FE, 1 },

    // FSCA FPUL, DRn - sine/cosine table lookup
    // TODO: the issue cycle count here might be wrong, I couldn't find that
    //       value for this instruction
    { "1111nnn011111101", FPU_HANDLER(fsca_fpu), false, SH4_GROUP_FE, 1 },

    // FSRRA FRn
    // 1111nnnn01111101
    // TODO: the issue cycle for this opcode might be wrong as well
    { "1111nnnn01111101", FPU_HANDLER(fsrra_fpu), false, SH4_GROUP_FE, 1 },

    { NULL }
};

static InstOpcode invalid_opcode = {
    "0000000000000000", &sh4_inst_invalid, false,
    (sh4_inst_group_t)0, 0, 0, 0
};

InstOpcode const *sh4_inst_lut[1 << 16];

void sh4_init_inst_lut() {
    unsigned inst;
    for (inst = 0; inst < (1 << 16); inst++)
        sh4_inst_lut[inst] = sh4_decode_inst_slow((inst_t)inst);
}

// used to initialize the sh4_inst_lut
static InstOpcode const* sh4_decode_inst_slow(inst_t inst) {
    InstOpcode const *op = opcode_list;

    while (op->fmt) {
        if ((op->mask & inst) == op->val) {
            return op;
        }
        op++;
    }

    return &invalid_opcode;
}

void sh4_compile_instructions(Sh4 *sh4) {
    InstOpcode *op = opcode_list;

    while (op->fmt) {
        sh4_compile_instruction(sh4, op);
        op++;
    }
}

void sh4_compile_instruction(Sh4 *sh4, struct InstOpcode *op) {
    char const *fmt = op->fmt;
    inst_t mask = 0, val = 0;

    if (strlen(fmt) != 16) {
        error_set_param_name("instruction opcode format");
        error_set_opcode_format(fmt);
        SH4_INST_RAISE_ERROR(sh4, ERROR_INVALID_PARAM);
    }

    for (int idx = 0; idx < 16; idx++) {
        val <<= 1;
        mask <<= 1;

        if (fmt[idx] == '1' || fmt[idx] == '0') {
            mask |= 1;
        }

        if (fmt[idx] == '1')
            val |= 1;
    }

    op->mask = mask;
    op->val = val;
}
