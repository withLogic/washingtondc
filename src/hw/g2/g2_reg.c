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

#include <string.h>
#include <stdio.h>

#include <stdint.h>

#include "types.h"
#include "mem_code.h"
#include "mem_areas.h"
#include "error.h"
#include "log.h"
#include "mmio.h"

#include "g2_reg.h"

#define N_G2_REGS (ADDR_G2_LAST - ADDR_G2_FIRST + 1)

DECL_MMIO_REGION(g2_reg_32, N_G2_REGS, ADDR_G2_FIRST, uint32_t)
DEF_MMIO_REGION(g2_reg_32, N_G2_REGS, ADDR_G2_FIRST, uint32_t)

static uint8_t reg_backing[N_G2_REGS];

uint8_t g2_reg_read_8(addr32_t addr) {
    error_set_length(1);
    error_set_address(addr);
    RAISE_ERROR(ERROR_UNIMPLEMENTED);
}

void g2_reg_write_8(addr32_t addr, uint8_t val) {
    error_set_length(1);
    error_set_address(addr);
    RAISE_ERROR(ERROR_UNIMPLEMENTED);
}

uint16_t g2_reg_read_16(addr32_t addr) {
    error_set_length(2);
    error_set_address(addr);
    RAISE_ERROR(ERROR_UNIMPLEMENTED);
}

void g2_reg_write_16(addr32_t addr, uint16_t val) {
    error_set_length(2);
    error_set_address(addr);
    RAISE_ERROR(ERROR_UNIMPLEMENTED);
}

uint32_t g2_reg_read_32(addr32_t addr) {
    return mmio_region_g2_reg_32_read(&mmio_region_g2_reg_32, addr);
}

void g2_reg_write_32(addr32_t addr, uint32_t val) {
    mmio_region_g2_reg_32_write(&mmio_region_g2_reg_32, addr, val);
}

float g2_reg_read_float(addr32_t addr) {
    uint32_t tmp = mmio_region_g2_reg_32_read(&mmio_region_g2_reg_32, addr);
    float ret;
    memcpy(&ret, &tmp, sizeof(ret));
    return ret;
}

void g2_reg_write_float(addr32_t addr, float val) {
    uint32_t tmp;
    memcpy(&tmp, &val, sizeof(tmp));
    mmio_region_g2_reg_32_write(&mmio_region_g2_reg_32, addr, tmp);
}

double g2_reg_read_double(addr32_t addr) {
    error_set_length(8);
    error_set_address(addr);
    RAISE_ERROR(ERROR_UNIMPLEMENTED);
}

void g2_reg_write_double(addr32_t addr, double val) {
    error_set_length(8);
    error_set_address(addr);
    RAISE_ERROR(ERROR_UNIMPLEMENTED);
}

static void sb_adst_reg_mmio_write(struct mmio_region_g2_reg_32 *region,
                                   unsigned idx, uint32_t val) {
    if (val) {
        error_set_feature("AICA DMA");
        RAISE_ERROR(ERROR_UNIMPLEMENTED);
    }
}

void g2_reg_init(void) {
    init_mmio_region_g2_reg_32(&mmio_region_g2_reg_32, (void*)reg_backing);

    mmio_region_g2_reg_32_init_cell(&mmio_region_g2_reg_32,
                                    "SB_ADSTAG", 0x5f7800,
                                    mmio_region_g2_reg_32_warn_read_handler,
                                    mmio_region_g2_reg_32_warn_write_handler);
    mmio_region_g2_reg_32_init_cell(&mmio_region_g2_reg_32,
                                    "SB_ADSTAR", 0x5f7804,
                                    mmio_region_g2_reg_32_warn_read_handler,
                                    mmio_region_g2_reg_32_warn_write_handler);
    mmio_region_g2_reg_32_init_cell(&mmio_region_g2_reg_32,
                                    "SB_ADLEN", 0x5f7808,
                                    mmio_region_g2_reg_32_warn_read_handler,
                                    mmio_region_g2_reg_32_warn_write_handler);
    mmio_region_g2_reg_32_init_cell(&mmio_region_g2_reg_32,
                                    "SB_ADDIR", 0x5f780c,
                                    mmio_region_g2_reg_32_warn_read_handler,
                                    mmio_region_g2_reg_32_warn_write_handler);
    mmio_region_g2_reg_32_init_cell(&mmio_region_g2_reg_32,
                                    "SB_ADTSEL", 0x5f7810,
                                    mmio_region_g2_reg_32_warn_read_handler,
                                    mmio_region_g2_reg_32_warn_write_handler);
    mmio_region_g2_reg_32_init_cell(&mmio_region_g2_reg_32,
                                    "SB_ADEN", 0x5f7814,
                                    mmio_region_g2_reg_32_warn_read_handler,
                                    mmio_region_g2_reg_32_warn_write_handler);
    mmio_region_g2_reg_32_init_cell(&mmio_region_g2_reg_32,
                                    "SB_ADST", 0x5f7818,
                                    mmio_region_g2_reg_32_warn_read_handler,
                                    sb_adst_reg_mmio_write);
    mmio_region_g2_reg_32_init_cell(&mmio_region_g2_reg_32,
                                    "SB_ADSUSP", 0x5f781c,
                                    mmio_region_g2_reg_32_warn_read_handler,
                                    mmio_region_g2_reg_32_warn_write_handler);
    mmio_region_g2_reg_32_init_cell(&mmio_region_g2_reg_32,
                                    "SB_E1STAG", 0x5f7820,
                                    mmio_region_g2_reg_32_warn_read_handler,
                                    mmio_region_g2_reg_32_warn_write_handler);
    mmio_region_g2_reg_32_init_cell(&mmio_region_g2_reg_32,
                                    "SB_E1STAR", 0x5f7824,
                                    mmio_region_g2_reg_32_warn_read_handler,
                                    mmio_region_g2_reg_32_warn_write_handler);
    mmio_region_g2_reg_32_init_cell(&mmio_region_g2_reg_32,
                                    "SB_E1LEN", 0x5f7828,
                                    mmio_region_g2_reg_32_warn_read_handler,
                                    mmio_region_g2_reg_32_warn_write_handler);
    mmio_region_g2_reg_32_init_cell(&mmio_region_g2_reg_32,
                                    "SB_E1DIR", 0x5f782c,
                                    mmio_region_g2_reg_32_warn_read_handler,
                                    mmio_region_g2_reg_32_warn_write_handler);
    mmio_region_g2_reg_32_init_cell(&mmio_region_g2_reg_32,
                                    "SB_E1TSEL", 0x5f7830,
                                    mmio_region_g2_reg_32_warn_read_handler,
                                    mmio_region_g2_reg_32_warn_write_handler);
    mmio_region_g2_reg_32_init_cell(&mmio_region_g2_reg_32,
                                    "SB_E1EN", 0x5f7834,
                                    mmio_region_g2_reg_32_warn_read_handler,
                                    mmio_region_g2_reg_32_warn_write_handler);
    mmio_region_g2_reg_32_init_cell(&mmio_region_g2_reg_32,
                                    "SB_E1ST", 0x5f7838,
                                    mmio_region_g2_reg_32_warn_read_handler,
                                    mmio_region_g2_reg_32_warn_write_handler);
    mmio_region_g2_reg_32_init_cell(&mmio_region_g2_reg_32,
                                    "SB_E1SUSP", 0x5f783c,
                                    mmio_region_g2_reg_32_warn_read_handler,
                                    mmio_region_g2_reg_32_warn_write_handler);
    mmio_region_g2_reg_32_init_cell(&mmio_region_g2_reg_32,
                                    "SB_E2STAG", 0x5f7840,
                                    mmio_region_g2_reg_32_warn_read_handler,
                                    mmio_region_g2_reg_32_warn_write_handler);
    mmio_region_g2_reg_32_init_cell(&mmio_region_g2_reg_32,
                                    "SB_E2STAR", 0x5f7844,
                                    mmio_region_g2_reg_32_warn_read_handler,
                                    mmio_region_g2_reg_32_warn_write_handler);
    mmio_region_g2_reg_32_init_cell(&mmio_region_g2_reg_32,
                                    "SB_E2LEN", 0x5f7848,
                                    mmio_region_g2_reg_32_warn_read_handler,
                                    mmio_region_g2_reg_32_warn_write_handler);
    mmio_region_g2_reg_32_init_cell(&mmio_region_g2_reg_32,
                                    "SB_E2DIR", 0x5f784c,
                                    mmio_region_g2_reg_32_warn_read_handler,
                                    mmio_region_g2_reg_32_warn_write_handler);
    mmio_region_g2_reg_32_init_cell(&mmio_region_g2_reg_32,
                                    "SB_E2TSEL", 0x5f7850,
                                    mmio_region_g2_reg_32_warn_read_handler,
                                    mmio_region_g2_reg_32_warn_write_handler);
    mmio_region_g2_reg_32_init_cell(&mmio_region_g2_reg_32,
                                    "SB_E2EN", 0x5f7854,
                                    mmio_region_g2_reg_32_warn_read_handler,
                                    mmio_region_g2_reg_32_warn_write_handler);
    mmio_region_g2_reg_32_init_cell(&mmio_region_g2_reg_32,
                                    "SB_E2ST", 0x5f7858,
                                    mmio_region_g2_reg_32_warn_read_handler,
                                    mmio_region_g2_reg_32_warn_write_handler);
    mmio_region_g2_reg_32_init_cell(&mmio_region_g2_reg_32,
                                    "SB_E2SUSP", 0x5f785c,
                                    mmio_region_g2_reg_32_warn_read_handler,
                                    mmio_region_g2_reg_32_warn_write_handler);
    mmio_region_g2_reg_32_init_cell(&mmio_region_g2_reg_32,
                                    "SB_DDSTAG", 0x5f7860,
                                    mmio_region_g2_reg_32_warn_read_handler,
                                    mmio_region_g2_reg_32_warn_write_handler);
    mmio_region_g2_reg_32_init_cell(&mmio_region_g2_reg_32,
                                    "SB_DDSTAR", 0x5f7864,
                                    mmio_region_g2_reg_32_warn_read_handler,
                                    mmio_region_g2_reg_32_warn_write_handler);
    mmio_region_g2_reg_32_init_cell(&mmio_region_g2_reg_32,
                                    "SB_DDLEN", 0x5f7868,
                                    mmio_region_g2_reg_32_warn_read_handler,
                                    mmio_region_g2_reg_32_warn_write_handler);
    mmio_region_g2_reg_32_init_cell(&mmio_region_g2_reg_32,
                                    "SB_DDIR", 0x5f786c,
                                    mmio_region_g2_reg_32_warn_read_handler,
                                    mmio_region_g2_reg_32_warn_write_handler);
    mmio_region_g2_reg_32_init_cell(&mmio_region_g2_reg_32,
                                    "SB_DDTSEL", 0x5f7870,
                                    mmio_region_g2_reg_32_warn_read_handler,
                                    mmio_region_g2_reg_32_warn_write_handler);
    mmio_region_g2_reg_32_init_cell(&mmio_region_g2_reg_32,
                                    "SB_DDEN", 0x5f7874,
                                    mmio_region_g2_reg_32_warn_read_handler,
                                    mmio_region_g2_reg_32_warn_write_handler);
    mmio_region_g2_reg_32_init_cell(&mmio_region_g2_reg_32,
                                    "SB_DDST", 0x5f7878,
                                    mmio_region_g2_reg_32_warn_read_handler,
                                    mmio_region_g2_reg_32_warn_write_handler);
    mmio_region_g2_reg_32_init_cell(&mmio_region_g2_reg_32,
                                    "SB_DDSUSP", 0x5f787c,
                                    mmio_region_g2_reg_32_warn_read_handler,
                                    mmio_region_g2_reg_32_warn_write_handler);

    /* some debugging bullshit, hopefully I never need these... */
    mmio_region_g2_reg_32_init_cell(&mmio_region_g2_reg_32,
                                    "SB_G2DSTO", 0x5f7890,
                                    mmio_region_g2_reg_32_warn_read_handler,
                                    mmio_region_g2_reg_32_warn_write_handler);
    mmio_region_g2_reg_32_init_cell(&mmio_region_g2_reg_32,
                                    "SB_G2TRTO", 0x5f7894,
                                    mmio_region_g2_reg_32_warn_read_handler,
                                    mmio_region_g2_reg_32_warn_write_handler);

    /* the modem, it will be a long time before I get around to this */
    mmio_region_g2_reg_32_init_cell(&mmio_region_g2_reg_32,
                                    "SB_G2MDMTO", 0x5f7898,
                                    mmio_region_g2_reg_32_warn_read_handler,
                                    mmio_region_g2_reg_32_warn_write_handler);
    mmio_region_g2_reg_32_init_cell(&mmio_region_g2_reg_32,
                                    "SB_G2MDMW", 0x5f789c,
                                    mmio_region_g2_reg_32_warn_read_handler,
                                    mmio_region_g2_reg_32_warn_write_handler);

    /* ??? */
    mmio_region_g2_reg_32_init_cell(&mmio_region_g2_reg_32,
                                    "UNKNOWN_G2_REG_0x5f78a0", 0x5f78a0,
                                    mmio_region_g2_reg_32_warn_read_handler,
                                    mmio_region_g2_reg_32_warn_write_handler);
    mmio_region_g2_reg_32_init_cell(&mmio_region_g2_reg_32,
                                    "UNKNOWN_G2_REG_0x5f78a4", 0x5f78a4,
                                    mmio_region_g2_reg_32_warn_read_handler,
                                    mmio_region_g2_reg_32_warn_write_handler);
    mmio_region_g2_reg_32_init_cell(&mmio_region_g2_reg_32,
                                    "UNKNOWN_G2_REG_0x5f78a8", 0x5f78a8,
                                    mmio_region_g2_reg_32_warn_read_handler,
                                    mmio_region_g2_reg_32_warn_write_handler);
    mmio_region_g2_reg_32_init_cell(&mmio_region_g2_reg_32,
                                    "UNKNOWN_G2_REG_0x5f78ac", 0x5f78ac,
                                    mmio_region_g2_reg_32_warn_read_handler,
                                    mmio_region_g2_reg_32_warn_write_handler);
    mmio_region_g2_reg_32_init_cell(&mmio_region_g2_reg_32,
                                    "UNKNOWN_G2_REG_0x5f78b0", 0x5f78b0,
                                    mmio_region_g2_reg_32_warn_read_handler,
                                    mmio_region_g2_reg_32_warn_write_handler);
    mmio_region_g2_reg_32_init_cell(&mmio_region_g2_reg_32,
                                    "UNKNOWN_G2_REG_0x5f78b4", 0x5f78b4,
                                    mmio_region_g2_reg_32_warn_read_handler,
                                    mmio_region_g2_reg_32_warn_write_handler);
    mmio_region_g2_reg_32_init_cell(&mmio_region_g2_reg_32,
                                    "UNKNOWN_G2_REG_0x5f78b8", 0x5f78b8,
                                    mmio_region_g2_reg_32_warn_read_handler,
                                    mmio_region_g2_reg_32_warn_write_handler);

    mmio_region_g2_reg_32_init_cell(&mmio_region_g2_reg_32,
                                    "SB_G2APRO", 0x5f78bc,
                                    mmio_region_g2_reg_32_warn_read_handler,
                                    mmio_region_g2_reg_32_warn_write_handler);
}

void g2_reg_cleanup(void) {
    cleanup_mmio_region_g2_reg_32(&mmio_region_g2_reg_32);
}
