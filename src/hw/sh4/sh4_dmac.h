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

#ifndef SH4_DMAC_H_
#define SH4_DMAC_H_

struct Sh4;
#include "types.h"

struct sh4_dmac {
    
};

int sh4_dmac_sar_reg_read_handler(Sh4 *sh4, void *buf,
                                  struct Sh4MemMappedReg const *reg_info);
int sh4_dmac_sar_reg_write_handler(Sh4 *sh4, void const *buf,
                                   struct Sh4MemMappedReg const *reg_info);
int sh4_dmac_dar_reg_read_handler(Sh4 *sh4, void *buf,
                                  struct Sh4MemMappedReg const *reg_info);
int sh4_dmac_dar_reg_write_handler(Sh4 *sh4, void const *buf,
                                   struct Sh4MemMappedReg const *reg_info);
int sh4_dmac_dmatcr_reg_read_handler(Sh4 *sh4, void *buf,
                                     struct Sh4MemMappedReg const *reg_info);
int sh4_dmac_dmatcr_reg_write_handler(Sh4 *sh4, void const *buf,
                                      struct Sh4MemMappedReg const *reg_info);
int sh4_dmac_chcr_reg_read_handler(Sh4 *sh4, void *buf,
                                   struct Sh4MemMappedReg const *reg_info);
int sh4_dmac_chcr_reg_write_handler(Sh4 *sh4, void const *buf,
                                    struct Sh4MemMappedReg const *reg_info);
int sh4_dmac_dmaor_reg_read_handler(Sh4 *sh4, void *buf,
                                    struct Sh4MemMappedReg const *reg_info);
int sh4_dmac_dmaor_reg_write_handler(Sh4 *sh4, void const *buf,
                                     struct Sh4MemMappedReg const *reg_info);

/*
 * perform a DMA transfer from some external device to memory.
 * This completes the transfer immediately instead of
 * modeling the cycle-steal/burst transfer characteristics.
 *
 * this function does not raise any interrupts.
 */
void sh4_dmac_transfer_to_mem(addr32_t transfer_dst, size_t unit_sz,
                              size_t n_units, void const *dat);

#endif