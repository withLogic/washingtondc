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

#ifndef CODE_BLOCK_X86_64_H_
#define CODE_BLOCK_X86_64_H_

#ifndef ENABLE_JIT_X86_64
#error this file should not be built when the x86_64 JIT backend is disabled
#endif

struct il_code_block;

struct code_block_x86_64 {
    /* void(*native)(void); */
    void *native;
    unsigned cycle_count;
    unsigned bytes_used;
};

void code_block_x86_64_init(struct code_block_x86_64 *blk);
void code_block_x86_64_cleanup(struct code_block_x86_64 *blk);

void code_block_x86_64_compile(struct code_block_x86_64 *out,
                               struct il_code_block const *il_blk);

#endif
