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

#ifndef MAPLE_H_
#define MAPLE_H_

#include <stdint.h>
#include <stdbool.h>

enum maple_cmd {
    // maplebus response codes
    MAPLE_RESP_NONE = -1,
    MAPLE_RESP_DEVINFO = 5,
    MAPLE_RESP_DATATRF = 8,

    // maplebus command codes
    MAPLE_CMD_DEVINFO = 1,
    MAPLE_CMD_GETCOND = 9
};

#define MAPLE_FRAME_OUTPUT_DATA_LEN 1024
#define MAPLE_FRAME_INPUT_DATA_LEN 1024

struct maple_frame {
    /* enum maple_cmd cmd; */
    unsigned port;
    unsigned ptrn;
    uint32_t recv_addr;

    bool last_frame;

    enum maple_cmd cmd;
    unsigned maple_addr;
    unsigned pack_len;

    unsigned input_len; // length of input in bytes (not words)
    uint8_t input_data[MAPLE_FRAME_INPUT_DATA_LEN];

    unsigned output_len; // length of input in bytes (not words)
    uint8_t output_data[MAPLE_FRAME_OUTPUT_DATA_LEN];
};

void maple_handle_frame(struct maple_frame *frame);

uint32_t maple_read_frame(struct maple_frame *frame_out, uint32_t addr);
void maple_write_frame_resp(struct maple_frame *frame, unsigned resp_code);

#define MAPLE_TRACE(msg, ...) LOG_DBG("MAPLE: "msg, ##__VA_ARGS__)

// don't call this directly, use the MAPLE_TRACE macro instead
void maple_do_trace(char const *msg, ...);

#define MAPLE_PORT_COUNT 4
#define MAPLE_UNIT_COUNT 6

unsigned maple_addr_pack(unsigned port, unsigned unit);
void maple_addr_unpack(unsigned addr, unsigned *port_out, unsigned *unit_out);

void maple_init(void);
void maple_cleanup(void);

#endif
