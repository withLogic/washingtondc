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

#include <stdlib.h>

#include "sh4_excp.h"
#include "sh4_mem.h"
#include "sh4_ocache.h"
#include "sh4.h"
#include "dreamcast.h"
#include "MemoryMap.h"
#include "mem_code.h"

#ifdef ENABLE_DEBUGGER
#include "debugger.h"
#endif

static inline enum VirtMemArea sh4_get_mem_area(addr32_t addr);

/*
 * TODO: need to adequately return control to the debugger when there's a memory
 * error and the debugger has its error-handler set up.  longjmp is the obvious
 * solution, but until all the codebase is out of C++ I don't want to risk that.
 */

#define SH4_WRITE_MEM_TMPL(type, postfix)                               \
    void sh4_write_mem_##postfix(Sh4 *sh4, type val, addr32_t addr) {   \
        enum VirtMemArea virt_area = sh4_get_mem_area(addr);            \
        switch (virt_area) {                                            \
        case SH4_AREA_P0:                                               \
        case SH4_AREA_P3:                                               \
            /*                                                          \
             * TODO: Check for MMUCR_AT_MASK in the MMUCR register and raise \
             * an error or do TLB lookups accordingly.                  \
             *                                                          \
             * currently it is impossible for this to be set because of the \
             * ERROR_UNIMPLEMENTED that gets raised if you set this bit in \
             * sh4_reg.c                                                \
             */                                                         \
                                                                        \
            /* handle the case where OCE is enabled and ORA is */       \
            /* enabled but we don't have Ocache available */            \
            if ((sh4->reg[SH4_REG_CCR] & SH4_CCR_OCE_MASK) &&           \
                (sh4->reg[SH4_REG_CCR] & SH4_CCR_ORA_MASK) &&           \
                sh4_ocache_in_ram_area(addr)) {                         \
                sh4_ocache_do_write_ora(sh4, &val, addr, sizeof(val));  \
                return;                                                 \
            }                                                           \
                                                                        \
            /* don't use the cache */                                   \
            /* INTENTIONAL FALLTHROUGH */                               \
        case SH4_AREA_P1:                                               \
        case SH4_AREA_P2:                                               \
            memory_map_write_##postfix(val, addr & 0x1fffffff);         \
            return;                                                     \
        case SH4_AREA_P4:                                               \
            if (sh4_do_write_p4(sh4, &val, addr, sizeof(val)) ==        \
                MEM_ACCESS_SUCCESS)                                     \
                return;                                                 \
            else                                                        \
                RAISE_ERROR(get_error_pending());                       \
        default:                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        error_set_wtf("this should not be possible");                   \
        RAISE_ERROR(ERROR_INTEGRITY);                                   \
        exit(1); /* never happens */                                    \
    }                                                                   \

SH4_WRITE_MEM_TMPL(uint8_t, 8)
SH4_WRITE_MEM_TMPL(uint16_t, 16)
SH4_WRITE_MEM_TMPL(uint32_t, 32)
SH4_WRITE_MEM_TMPL(float, float)
SH4_WRITE_MEM_TMPL(double, double)

#define SH4_READ_MEM_TMPL(type, postfix)                                \
    type sh4_read_mem_##postfix(Sh4 *sh4, addr32_t addr) {              \
        type tmp_val;                                                   \
        enum VirtMemArea virt_area = sh4_get_mem_area(addr);            \
        switch (virt_area) {                                            \
        case SH4_AREA_P0:                                               \
        case SH4_AREA_P3:                                               \
            /*                                                          \
             * TODO: Check for MMUCR_AT_MASK in the MMUCR register and raise \
             * an error or do TLB lookups accordingly.                  \
             *                                                          \
             * currently it is impossible for this to be set because of the \
             * ERROR_UNIMPLEMENTED that gets raised if you set this bit in \
             * sh4_reg.c                                                \
             */                                                         \
                                                                        \
            /* handle the case where OCE is enabled and ORA is */       \
            /* enabled but we don't have Ocache available */            \
            if ((sh4->reg[SH4_REG_CCR] & SH4_CCR_OCE_MASK) &&           \
                (sh4->reg[SH4_REG_CCR] & SH4_CCR_ORA_MASK) &&           \
                sh4_ocache_in_ram_area(addr)) {                         \
                sh4_ocache_do_read_ora(sh4, &tmp_val,                   \
                                       addr, sizeof(tmp_val));          \
                return tmp_val;                                         \
            }                                                           \
                                                                        \
            /* don't use the cache */                                   \
            /* INTENTIONAL FALLTHROUGH */                               \
        case SH4_AREA_P1:                                               \
        case SH4_AREA_P2:                                               \
            return memory_map_read_##postfix(addr & 0x1fffffff);        \
        case SH4_AREA_P4:                                               \
            if (sh4_do_read_p4(sh4, &tmp_val, addr, sizeof(tmp_val)) == MEM_ACCESS_SUCCESS) \
                return tmp_val;                                         \
            RAISE_ERROR(get_error_pending());                           \
        default:                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        /* TODO: memory access exception ? */                           \
        RAISE_ERROR(ERROR_UNIMPLEMENTED);                               \
    }

SH4_READ_MEM_TMPL(uint8_t, 8);
SH4_READ_MEM_TMPL(uint16_t, 16);
SH4_READ_MEM_TMPL(uint32_t, 32);
SH4_READ_MEM_TMPL(float, float);
SH4_READ_MEM_TMPL(double, double);

int sh4_do_read_p4(Sh4 *sh4, void *dat, addr32_t addr, unsigned len) {
    if ((addr & SH4_SQ_AREA_MASK) == SH4_SQ_AREA_VAL)
        return sh4_sq_read(sh4, dat, addr, len);

    if (addr >= SH4_P4_REGSTART && addr < SH4_P4_REGEND) {
        return sh4_read_mem_mapped_reg(sh4, dat, addr, len);
    }

    if (addr >= SH4_OC_ADDR_ARRAY_FIRST && addr <= SH4_OC_ADDR_ARRAY_LAST) {
        sh4_ocache_read_addr_array(sh4, dat, addr, len);
        return MEM_ACCESS_SUCCESS;
    }

    error_set_address(addr);
    error_set_feature("reading from part of the P4 memory region");
    PENDING_ERROR(ERROR_UNIMPLEMENTED);
    return MEM_ACCESS_FAILURE;
}

int sh4_do_write_p4(Sh4 *sh4, void const *dat, addr32_t addr, unsigned len) {
    if ((addr & SH4_SQ_AREA_MASK) == SH4_SQ_AREA_VAL)
        return sh4_sq_write(sh4, dat, addr, len);

    if (addr >= SH4_P4_REGSTART && addr < SH4_P4_REGEND) {
        return sh4_write_mem_mapped_reg(sh4, dat, addr, len);
    }

    if (addr >= SH4_OC_ADDR_ARRAY_FIRST && addr <= SH4_OC_ADDR_ARRAY_LAST) {
        sh4_ocache_write_addr_array(sh4, dat, addr, len);
        return MEM_ACCESS_SUCCESS;
    }

    error_set_address(addr);
    error_set_feature("writing to part of the P4 memory region");
    PENDING_ERROR(ERROR_UNIMPLEMENTED);
    return MEM_ACCESS_FAILURE;
}


static inline enum VirtMemArea sh4_get_mem_area(addr32_t addr) {
    /*
     * XXX I tried replacing this block of if statements with a lookup table,
     * but somehow it turned out to be slower that way.  This is possibly
     * because the lookup-table was not in the cache and had to be fetched from
     * memory.
     *
     * If you ever want to look into this again, the trick is to use the upper
     * four bits as the index into the lookup table (P0 will be 0-7,
     * P1 will be 8-9, etc.)
     */
    if (addr <= SH4_AREA_P0_LAST)
        return SH4_AREA_P0;
    if (addr <= SH4_AREA_P1_LAST)
        return SH4_AREA_P1;
    if (addr <= SH4_AREA_P2_LAST)
        return SH4_AREA_P2;
    if (addr <= SH4_AREA_P3_LAST)
        return SH4_AREA_P3;
    return SH4_AREA_P4;
}
