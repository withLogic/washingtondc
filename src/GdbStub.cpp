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

#include <iostream>
#include <iomanip>

#include "hw/sh4/sh4.hpp"
#include "Dreamcast.hpp"

#include "gdb_stub.h"

// uncomment this to log all traffic in/out of the debugger to stdout
// #define GDBSTUB_VERBOSE

static void gdb_on_break(void *arg);
static void gdb_on_read_watchpoint(addr32_t addr, void *arg);
static void gdb_on_write_watchpoint(addr32_t addr, void *arg);
static void gdb_on_softbreak(inst_t inst, addr32_t addr, void *arg);
static int decode_hex(char ch);

static void craft_packet(struct string *out, struct string const *in);
static void gdb_serialize_regs(struct gdb_stub *stub, struct string *out);
static void deserialize_regs(struct string const *input_str,
                             reg32_t regs[N_REGS]);
static void serialize_data(struct string *out, void const *buf,
                           unsigned buf_len);
static size_t deserialize_data(struct string const *input_str,
                               void *out, size_t max_sz);
static int decode_hex(char ch);
static void write_start(struct gdb_stub *stub);
static void extract_packet(struct string *out, struct string const *packet_in);
static int set_reg(reg32_t reg_file[SH4_REGISTER_COUNT], Sh4::FpuReg *fpu,
                   unsigned reg_no, reg32_t reg_val, bool bank);
static void handle_packet(struct gdb_stub *stub, struct string *pkt);
static void transmit_pkt(struct gdb_stub *stub, struct string const *pkt);
static bool next_packet(struct gdb_stub *stub, struct string *pkt);
// static void errror_handler(int error_tp, void *argptr);
static void expect_mem_access_error(struct gdb_stub *stub, bool should);

static void handle_c_packet(struct gdb_stub *stub, struct string *out,
                            struct string *dat);
static void handle_q_packet(struct gdb_stub *stub, struct string *out,
                            struct string const *dat);
static void handle_g_packet(struct gdb_stub *stub, struct string *out,
                            struct string const *dat);
static void handle_m_packet(struct gdb_stub *stub, struct string *out,
                            struct string const *dat);
static void handle_M_packet(struct gdb_stub *stub, struct string *out,
                            struct string const *dat);
static void handle_s_packet(struct gdb_stub *stub, struct string *out,
                            struct string const *dat);
static void handle_G_packet(struct gdb_stub *stub, struct string *out,
                            struct string const *dat);
static void handle_P_packet(struct gdb_stub *stub, struct string *out,
                            struct string const *dat);
static void handle_D_packet(struct gdb_stub *stub, struct string *out,
                            struct string const *dat);
static void handle_K_packet(struct gdb_stub *stub, struct string *out,
                            struct string const *dat);
static void handle_Z_packet(struct gdb_stub *stub, struct string *out,
                            struct string const *dat);
static void handle_z_packet(struct gdb_stub *stub, struct string *out,
                            struct string const *dat);

static void read_mem_4(struct gdb_stub *stub, struct string *out,
                       addr32_t addr, unsigned len);
static void read_mem_2(struct gdb_stub *stub, struct string *out,
                       addr32_t addr, unsigned len);
static void read_mem_1(struct gdb_stub *stub, struct string *out,
                       addr32_t addr, unsigned len);

static void
listener_cb(struct evconnlistener *listener,
            evutil_socket_t fd, struct sockaddr *saddr,
            int socklen, void *arg);
static void handle_events(struct bufferevent *bev, short events, void *arg);
static void handle_read(struct bufferevent *bev, void *arg);
static void handle_write(struct bufferevent *bev, void *arg);

static size_t deserialize_data(struct string const *input_str,
                               void *out, size_t max_sz) {
    uint8_t *out8 = (uint8_t*)out;
    size_t bytes_written = 0;
    char ch;
    // char const eof = std::char_traits<char>::eof();
    char const *input = string_get(input_str);

    while ((ch = *input++)/* && (ch != eof)*/) {
        if (bytes_written >= max_sz)
            return max_sz;

        *out8 = uint8_t(decode_hex(ch));
        bytes_written++;

        if ((ch = *input++)) {
            *out8 <<= 4;
            *out8 |= uint8_t(decode_hex(ch));
        } else {
            break;
        }

        out8++;
    }

    return bytes_written;
}

extern "C" void gdb_init(struct gdb_stub *stub, struct debugger *dbg) {
    string_init(&stub->unack_packet);
    string_init(&stub->input_packet);

    stub->frontend_supports_swbreak = false;
    stub->is_writing = false;
    stub->should_expect_mem_access_error =false;
    stub->listener = NULL;
    stub->is_listening = false;
    stub->bev = NULL;
    stub->dbg = dbg;

    stub->output_buffer = evbuffer_new();
    if (!stub->output_buffer)
        RAISE_ERROR(ERROR_FAILED_ALLOC);

    dbg->frontend.step = NULL;
    dbg->frontend.attach = gdb_attach;
    dbg->frontend.on_break = gdb_on_break;
    dbg->frontend.on_read_watchpoint = gdb_on_read_watchpoint;
    dbg->frontend.on_write_watchpoint = gdb_on_write_watchpoint;
    dbg->frontend.on_softbreak = gdb_on_softbreak;
    dbg->frontend.arg = stub;
}

extern "C" void gdb_cleanup(struct gdb_stub *stub) {
    // TODO: cleanup

    string_cleanup(&stub->input_packet);
    string_cleanup(&stub->unack_packet);
}

extern "C" void gdb_attach(void *argptr) {
    struct gdb_stub *stub = (struct gdb_stub*)argptr;

    std::cout << "Awaiting remote GDB connection on port " << GDB_PORT_NO <<
        "..." << std::endl;

    struct sockaddr_in sin;
    memset(&sin, 0, sizeof(sin));
    sin.sin_family = AF_INET;
    sin.sin_port = htons(GDB_PORT_NO);
    unsigned event_flags = LEV_OPT_REUSEABLE | LEV_OPT_CLOSE_ON_FREE;
    stub->listener = evconnlistener_new_bind(dc_event_base, listener_cb, stub,
                                             event_flags, -1,
                                             (struct sockaddr*)&sin, sizeof(sin));

    if (!stub->listener)
        RAISE_ERROR(ERROR_FAILED_ALLOC);

    // the listener_cb will set is_listening = false when we have a connection
    stub->is_listening = true;
    do {
        std::cout << "still waiting..." << std::endl;
        if (event_base_loop(dc_event_base, EVLOOP_ONCE) != 0)
            exit(4);
    } while (stub->is_listening);

    std::cout << "Connection established." << std::endl;
}

static void gdb_on_break(void *arg) {
    struct gdb_stub *stub = (struct gdb_stub*)arg;

    struct string resp, pkt;
    string_init(&pkt);
    string_init_txt(&resp, "S05");

    craft_packet(&pkt, &resp);
    transmit_pkt(stub, &pkt);

    string_cleanup(&resp);
    string_cleanup(&pkt);
}

static void gdb_on_softbreak(inst_t inst, addr32_t addr, void *arg) {
    struct gdb_stub *stub = (struct gdb_stub*)arg;

    struct string resp, pkt;
    string_init(&pkt);
    string_init(&resp);

    if (stub->frontend_supports_swbreak) {
        string_set(&resp, "T05swbreak:");
        string_append_hex32(&resp, addr);
        string_append_char(&resp, ';');
    } else {
        string_set(&resp, "T05swbreak:");
    }

    stub->dbg->cur_state = DEBUG_STATE_BREAK;

    craft_packet(&pkt, &resp);
    transmit_pkt(stub, &pkt);
}

static void gdb_on_read_watchpoint(addr32_t addr, void *arg) {
    struct gdb_stub *stub = (struct gdb_stub*)arg;

    struct string resp, pkt;
    string_init(&pkt);
    string_init_txt(&resp, "S05");

    craft_packet(&pkt, &resp);
    transmit_pkt(stub, &pkt);
}

static void gdb_on_write_watchpoint(addr32_t addr, void *arg) {
    struct gdb_stub *stub = (struct gdb_stub*)arg;

    struct string resp, pkt;
    string_init(&pkt);
    string_init_txt(&resp, "S05");

    craft_packet(&pkt, &resp);
    transmit_pkt(stub, &pkt);
}

static void gdb_serialize_regs(struct gdb_stub *stub, struct string *out) {
    Sh4 *cpu = dreamcast_get_cpu();
    reg32_t reg_file[SH4_REGISTER_COUNT];
    sh4_get_regs(cpu, reg_file);
    Sh4::FpuReg fpu_reg = sh4_get_fpu(cpu);
    reg32_t regs[N_REGS] = { 0 };

    // general-purpose registers
    for (int i = 0; i < 16; i++) {
        if (i < 8) {
            if (reg_file[SH4_REG_SR] & SH4_SR_RB_MASK)
                regs[R0 + i] = reg_file[SH4_REG_R0_BANK1 + i];
            else
                regs[R0 + i] = reg_file[SH4_REG_R0_BANK0 + i];
        }else {
            regs[R0 + i] = reg_file[SH4_REG_R8 + (i - 8)];
        }
    }

    // banked registers
    for (int i = 0; i < 8; i++) {
        regs[R0B0 + i] = reg_file[SH4_REG_R0_BANK0 + i];
        regs[R0B1 + i] = reg_file[SH4_REG_R0_BANK1 + i];
    }

    // TODO: floating point registers

    // system/control registers
    regs[PC] = reg_file[SH4_REG_PC];
    regs[PR] = reg_file[SH4_REG_PR];
    regs[GBR] = reg_file[SH4_REG_GBR];
    regs[VBR] = reg_file[SH4_REG_VBR];
    regs[MACH] = reg_file[SH4_REG_MACH];
    regs[MACL] = reg_file[SH4_REG_MACL];
    regs[SR] = reg_file[SH4_REG_SR];
    regs[SSR] = reg_file[SH4_REG_SSR];
    regs[SPC] = reg_file[SH4_REG_SPC];

    // FPU system/control registers
    regs[FPUL] = fpu_reg.fpul;
    regs[FPSCR] = fpu_reg.fpscr;

    serialize_data(out, regs, sizeof(regs));
}

static void deserialize_regs(struct string const *input_str, reg32_t regs[N_REGS]) {
    size_t sz_expect = N_REGS * sizeof(*regs);

    size_t sz_actual = deserialize_data(input_str, regs, sz_expect);

    if (sz_expect != sz_actual) {
        // TODO: better error messages
        std::cout << "sz_expect is " << sz_expect << ", sz_actual is " <<
            sz_actual << std::endl;
        RAISE_ERROR(ERROR_INTEGRITY);
    }
}

static void serialize_data(struct string *out, void const *buf,
                           unsigned buf_len) {
    uint8_t const *buf8 = (uint8_t const*)buf;
    static const char hex_tbl[16] = {
        '0', '1', '2', '3',
        '4', '5', '6', '7',
        '8', '9', 'a', 'b',
        'c', 'd', 'e', 'f'
    };

    for (unsigned i = 0; i < buf_len; i++) {
        string_append_char(out, hex_tbl[(*buf8) >> 4]);
        string_append_char(out, hex_tbl[(*buf8) & 0xf]);
        buf8++;
    }
}

static int decode_hex(char ch)
{
    if ((ch >= 'a') && (ch <= 'f'))
        return ch - 'a' + 10;
    if ((ch >= '0') && (ch <= '9'))
        return ch - '0';
    if ((ch >= 'A') && (ch <= 'F'))
        return ch - 'A' + 10;
    return -1;
}

static void transmit(struct gdb_stub *stub, struct string const *data) {
    int len = string_length(data);
    if (len > 0) {
        char const *data_c_str = string_get(data);

        if (evbuffer_add(stub->output_buffer, data_c_str,
                         sizeof(char) * len) < 0) {
            RAISE_ERROR(ERROR_FAILED_ALLOC);
        }

        write_start(stub);
    }
}

static void write_start(struct gdb_stub *stub) {
    if (!evbuffer_get_length(stub->output_buffer))
        stub->is_writing = false;

    if (stub->is_writing || !evbuffer_get_length(stub->output_buffer))
        return;

    bufferevent_write_buffer(stub->bev, stub->output_buffer);
}

static void craft_packet(struct string *out, struct string const *in) {
    uint8_t csum = 0;
    static const char hex_tbl[16] = {
        '0', '1', '2', '3',
        '4', '5', '6', '7',
        '8', '9', 'a', 'b',
        'c', 'd', 'e', 'f'
    };

    char const *str_ptr = string_get(in);
    while (*str_ptr) {
        csum += (uint8_t)*str_ptr;
        str_ptr++;
    }

    string_append(out, "$");
    string_append(out, string_get(in));
    string_append(out, "#");
    string_append_char(out, hex_tbl[csum >> 4]);
    string_append_char(out, hex_tbl[csum & 0xf]);
}

static void extract_packet(struct string *out, struct string const *packet_in) {
    int dollar_idx = string_find_first_of(packet_in, "$");
    int pound_idx = string_find_last_of(packet_in, "#");

    if (dollar_idx < 0 || pound_idx < 0)
        return;

    string_substr(out, packet_in, dollar_idx + 1, pound_idx - 1);
}

static int set_reg(reg32_t reg_file[SH4_REGISTER_COUNT], Sh4::FpuReg *fpu,
                   unsigned reg_no, reg32_t reg_val, bool bank) {
    // there is some ambiguity over whether register banking should be based off
    // of the old sr or the new sr.  For now, it's based off of the old sr.

    // TODO: floating point registers
    if (reg_no >= R0 && reg_no <= R15) {
        unsigned idx = reg_no - R0;

        if (idx < 8) {
            if (bank)
                reg_file[SH4_REG_R0_BANK1 + idx] = reg_val;
            else
                reg_file[SH4_REG_R0_BANK0 + idx] = reg_val;
        } else {
            reg_file[SH4_REG_R8 + (idx + 8)] = reg_val;
        }
    } else if (reg_no >= R0B0 && reg_no <= R7B0) {
        reg_file[reg_no - R0B0 + SH4_REG_R0_BANK0] = reg_val;
    } else if (reg_no >= R0B1 && reg_no <= R7B1) {
        reg_file[reg_no - R0B1 + SH4_REG_R0_BANK1] = reg_val;
    } else if (reg_no == PC) {
        reg_file[SH4_REG_PC] =reg_val;
    } else if (reg_no == PR) {
        reg_file[SH4_REG_PR] = reg_val;
    } else if (reg_no == GBR) {
        reg_file[SH4_REG_GBR] = reg_val;
    } else if (reg_no == VBR) {
        reg_file[SH4_REG_VBR] = reg_val;
    } else if (reg_no == MACH) {
        reg_file[SH4_REG_MACH] = reg_val;
    } else if (reg_no == MACL) {
        reg_file[SH4_REG_MACL] = reg_val;
    } else if (reg_no == SR) {
        reg_file[SH4_REG_SR] = reg_val;
    } else if (reg_no == SSR) {
        reg_file[SH4_REG_SSR] = reg_val;
    } else if (reg_no == SPC) {
        reg_file[SH4_REG_SPC] = reg_val;
    } else if (reg_no == FPUL) {
        fpu->fpul = reg_val;
    } else if (reg_no == FPSCR) {
        fpu->fpscr = reg_val;
    } else {
#ifdef GDBSTUB_VERBOSE
        std::cout << "WARNING: GdbStub unable to set value of register " <<
            std::hex << reg_no << " to " << reg_val << std::endl;
#endif
        return 1;
    }

    return 0;
}

static void err_str(struct string *err_out, unsigned err_val) {
    static char const hex_chars[] = {
        '0', '1', '2', '3',
        '4', '5', '6', '7',
        '8', '9', 'a', 'b',
        'c', 'd', 'e', 'f'
    };

    string_append_char(err_out, 'E');

    // dont print more than 2 digits
    err_val &= 0xff;
    string_append_char(err_out, hex_chars[err_val >> 4]);
    err_val &= 0x0f;
    string_append_char(err_out, hex_chars[err_val]);
}

static void handle_c_packet(struct gdb_stub *stub, struct string *out,
                            struct string *dat) {
    stub->dbg->cur_state = DEBUG_STATE_NORM;
}

static void handle_q_packet(struct gdb_stub *stub, struct string *out,
                            struct string const *dat_orig) {
    struct string dat, tok;
    string_init(&dat);
    string_init(&tok);
    string_copy(&dat, dat_orig);

    if (string_eq_n(&dat, "qSupported", 10)) {
        int semicolon_idx = string_find_first_of(&dat, ";");

        if (semicolon_idx == -1)
            goto cleanup;

        struct string tmp;
        string_init(&tmp);
        string_substr(&tmp, &dat, semicolon_idx + 1, string_length(&dat) - 1);
        string_cleanup(&dat);
        memcpy(&dat, &tmp, sizeof(dat));

        struct string_curs curs;
        string_tok_begin(&curs);
        while (string_tok_next(&tok, &curs, string_get(&dat), ";")) {
            bool supported = false;

            int plus_or_minus_idx = string_find_last_of(&tok, "+-");

            /*
             * ignore all the settings that try to set variables,
             * we're really only here for swbreak.
             */
            if (plus_or_minus_idx >= 0) {
                if (string_get(&tok)[plus_or_minus_idx] == '+')
                    supported = true;

                struct string tmp;
                string_init(&tmp);
                string_substr(&tmp, &tok, 0, plus_or_minus_idx - 1);
                string_cleanup(&tok);
                memcpy(&tok, &tmp, sizeof(tok));
            }

            if (strcmp(string_get(&tok), "swbreak") == 0) {
                if (supported) {
                    stub->frontend_supports_swbreak = true;
                    string_append(out, "swbreak+;");
                } else {
                    string_append(out, "swbreak-;");
                }
            } else {
                string_append(out, string_get(&tok));
                string_append(out, "-;");
            }
        }

        goto cleanup;
    }

cleanup:
    string_cleanup(&tok);
    string_cleanup(&dat);
}

static void handle_g_packet(struct gdb_stub *stub, struct string *out,
                            struct string const *dat) {
    gdb_serialize_regs(stub, out);
}

static void handle_m_packet(struct gdb_stub *stub, struct string *out,
                            struct string const *dat) {
    int addr_idx = string_find_last_of(dat, "m");
    int comma_idx = string_find_last_of(dat, ",");
    int len_idx = comma_idx;

    if (addr_idx < 0 || comma_idx < 0 || len_idx < 0) {
        err_str(out, EINVAL);
        return;
    }

    addr_idx++;
    len_idx++;

    struct string len_str, addr_str;
    string_init(&len_str);
    string_init(&addr_str);

    string_substr(&len_str, dat, len_idx, string_length(dat) - 1);
    string_substr(&addr_str, dat, addr_idx, comma_idx - 1);

    uint32_t len = string_read_hex32(&len_str, 0);
    uint32_t addr = string_read_hex32(&addr_str, 0);

    string_cleanup(&addr_str);
    string_cleanup(&len_str);

    if (len % 4 == 0)
        return read_mem_4(stub, out, addr, len);
    else if (len % 2 == 0)
        return read_mem_2(stub, out, addr, len);
    else
        return read_mem_1(stub, out, addr, len);
}

static void read_mem_4(struct gdb_stub *stub, struct string *out,
                       addr32_t addr, unsigned len) {
    expect_mem_access_error(stub, true);

    while (len) {
        uint32_t val;

        try {
            sh4_read_mem(dreamcast_get_cpu(), &val, addr, sizeof(val));
            addr += 4;
        } catch (BaseException& exc) {
            expect_mem_access_error(stub, false);
            err_str(out, EINVAL);
            return;
        }

        if (stub->mem_access_error) {
            expect_mem_access_error(stub, false);
            err_str(out, EINVAL);
            return;
        }

        serialize_data(out, &val, sizeof(val));
        len -= 4;
    }

    expect_mem_access_error(stub, false);
}

static void read_mem_2(struct gdb_stub *stub, struct string *out,
                       addr32_t addr, unsigned len) {
    expect_mem_access_error(stub, true);

    while (len) {
        uint16_t val;

        try {
            sh4_read_mem(dreamcast_get_cpu(), &val, addr, sizeof(val));
            addr += 2;
        } catch (BaseException& exc) {
            expect_mem_access_error(stub, false);
            err_str(out, EINVAL);
            return;
        }

        if (stub->mem_access_error) {
            expect_mem_access_error(stub, false);
            err_str(out, EINVAL);
            return;
        }

        serialize_data(out, &val, sizeof(val));
        len -= 2;
    }

    expect_mem_access_error(stub, false);
}

static void read_mem_1(struct gdb_stub *stub, struct string *out,
                       addr32_t addr, unsigned len) {
    expect_mem_access_error(stub, true);

    while (len--) {
        uint8_t val;

        try {
            sh4_read_mem(dreamcast_get_cpu(), &val, addr, sizeof(val));
            addr++;
        } catch (BaseException& exc) {
            expect_mem_access_error(stub, false);
            err_str(out, EINVAL);
            return;
        }

        if (stub->mem_access_error) {
            expect_mem_access_error(stub, false);
            err_str(out, EINVAL);
            return;
        }

        serialize_data(out, &val, sizeof(val));
    }

    expect_mem_access_error(stub, false);
}

/*
 * TODO: bounds checking (not that I expect there to be any hackers going in
 * through the debugger of all places)
 */
static void handle_M_packet(struct gdb_stub *stub, struct string *out,
                            struct string const *dat) {
    int addr_idx = string_find_last_of(dat, "M");
    int comma_idx = string_find_last_of(dat, ",");
    int colon_idx = string_find_last_of(dat, ":");

    if (addr_idx < 0 || comma_idx < 0 || colon_idx < 0) {
        err_str(out, EINVAL);
        return;
    }

    int len_idx = comma_idx + 1;
    int dat_idx = colon_idx + 1;
    addr_idx++;

    struct string addr_substr;
    struct string len_substr;
    string_init(&addr_substr);
    string_init(&len_substr);
    string_substr(&addr_substr, dat, addr_idx, comma_idx - 1);
    string_substr(&len_substr, dat, len_idx, colon_idx - 1);

    uint32_t addr = string_read_hex32(&addr_substr, 0);
    uint32_t len = string_read_hex32(&len_substr, 0);

    string_cleanup(&addr_substr);
    string_cleanup(&len_substr);

    try {
        struct string new_dat;
        string_init(&new_dat);
        string_substr(&new_dat, dat, dat_idx, string_length(dat) - 1);
        if (len < 1024) {
            uint8_t *buf = new uint8_t[len];
            try {
                expect_mem_access_error(stub, true);
                deserialize_data(&new_dat, buf, len);
                sh4_write_mem(dreamcast_get_cpu(), buf, addr, len);

                if (stub->mem_access_error) {
                    expect_mem_access_error(stub, false);
                    err_str(out, EINVAL);
                    string_cleanup(&new_dat);
                    return;
                }

            } catch (BaseException& exc) {
                expect_mem_access_error(stub, false);
                string_cleanup(&new_dat);
                delete[] buf;
                throw;
            }
            delete[] buf;
        } else {
            error_set_length(len);
            RAISE_ERROR(ERROR_INVALID_PARAM);
        }
        string_cleanup(&new_dat);
    } catch (BaseException& exc) {
        expect_mem_access_error(stub, false);
        std::cerr << boost::diagnostic_information(exc);
        err_str(out, EINVAL);
        return;
    }

    expect_mem_access_error(stub, false);
    string_set(out, "OK");
}

static void handle_s_packet(struct gdb_stub *stub, struct string *out,
                            struct string const *dat) {
    stub->dbg->cur_state = DEBUG_STATE_PRE_STEP;
}

static void handle_G_packet(struct gdb_stub *stub, struct string *out,
                            struct string const *dat) {
    reg32_t regs[N_REGS];

    struct string tmp;
    string_init(&tmp);
    string_substr(&tmp, dat, 1, string_length(dat) - 1);
    deserialize_regs(&tmp, regs);
    string_cleanup(&tmp);

    reg32_t new_regs[SH4_REGISTER_COUNT];
    sh4_get_regs(dreamcast_get_cpu(), new_regs);
    Sh4::FpuReg new_fpu = sh4_get_fpu(dreamcast_get_cpu());
    bool bank = new_regs[SH4_REG_SR] & SH4_SR_RB_MASK;

    for (unsigned reg_no = 0; reg_no < N_REGS; reg_no++)
        set_reg(new_regs, &new_fpu, reg_no, regs[reg_no], bank);

    string_set(out, "OK");
}

static void handle_P_packet(struct gdb_stub *stub, struct string *out,
                            struct string const *dat) {
    Sh4 *cpu;
    reg32_t regs[SH4_REGISTER_COUNT];
    Sh4::FpuReg fpu;
    int equals_idx = string_find_first_of(dat, "=");

    if ((equals_idx < 0) || ((size_t)equals_idx >= (string_length(dat) - 1))) {
#ifdef GDBSTUB_VERBOSE
        std::cout << "WARNING: malformed P packet in gdbstub \"" << dat <<
            "\"" << std::endl;
#endif

        string_set(out, "E16");
        return;
    }

    struct string reg_no_str, reg_val_str;
    string_init(&reg_no_str);
    string_init(&reg_val_str);

    string_substr(&reg_no_str, dat, 1, equals_idx - 1);
    string_substr(&reg_val_str, dat, equals_idx + 1, string_length(dat) - 1);

    unsigned reg_no = 0;
    reg32_t reg_val = 0;
    deserialize_data(&reg_no_str, &reg_no, sizeof(reg_no));
    deserialize_data(&reg_val_str, &reg_val, sizeof(reg_val));

    if (reg_no >= N_REGS) {
#ifdef GDBSTUB_VERBOSE
        std::cout << "ERROR: unable to write to register number " <<
            std::hex << reg_no << std::endl;
#endif
        string_set(out, "E16");
        goto cleanup;
    }

    cpu = dreamcast_get_cpu();
    sh4_get_regs(cpu, regs);
    fpu = sh4_get_fpu(cpu);
    set_reg(regs, &fpu, reg_no, reg_val,
            bool(regs[SH4_REG_SR] & SH4_SR_RB_MASK));
    sh4_set_regs(cpu, regs);

    string_set(out, "OK");

cleanup:
    string_cleanup(&reg_val_str);
    string_cleanup(&reg_no_str);
}

static void handle_D_packet(struct gdb_stub *stub, struct string *out,
                            struct string const *dat) {
    stub->dbg->cur_state = DEBUG_STATE_NORM;

    debug_on_detach(stub->dbg);

    string_set(out, "OK");
}

static void handle_K_packet(struct gdb_stub *stub, struct string *out,
                            struct string const *dat) {
    dreamcast_kill();

    string_set(out, "OK");
}

static void handle_Z_packet(struct gdb_stub *stub, struct string *out,
                            struct string const *dat) {
    char const *dat_c_str = string_get(dat);
    struct string dat_local;
    string_init(&dat_local);
    string_copy(&dat_local, dat);

    if (dat_c_str[1] == '1') {
        //hardware breakpoint

        if (string_find_first_of(&dat_local, ":") >= 0)
            goto cleanup; // we don't support conditions

        int first_comma_idx = string_find_first_of(&dat_local, ",");
        if ((first_comma_idx < 0) ||
            (first_comma_idx == (int)string_length(&dat_local) - 1)) {
            goto cleanup; // something is wrong and/or unexpected
        }

        int last_comma_idx = string_find_last_of(&dat_local, ",");
        if (last_comma_idx < 0)
            goto cleanup; // something is wrong and/or unexpected

        struct string tmp_str;
        string_init(&tmp_str);
        string_substr(&tmp_str, &dat_local,
                      first_comma_idx + 1, last_comma_idx - 1);
        string_cleanup(&dat_local);
        memcpy(&dat_local, &tmp_str, sizeof(dat_local));

        addr32_t break_addr = string_read_hex32(&dat_local, 0);

        int err_code = debug_add_break(stub->dbg, break_addr);

        if (err_code == 0)
            string_set(out, "OK");
        else
            err_str(out, err_code);

        goto cleanup;
    } else if (dat_c_str[1] == '2') {
        // write watchpoint

        if (string_find_first_of(&dat_local, ":") >= 0)
            goto cleanup; // we don't support conditions

        int first_comma_idx = string_find_first_of(&dat_local, ",");
        if ((first_comma_idx < 0) ||
            (first_comma_idx == (int)string_length(&dat_local) - 1)) {
            goto cleanup; // something is wrong and/or unexpected
        }

        int last_comma_idx = string_find_last_of(&dat_local, ",");
        if (last_comma_idx < 0)
            goto cleanup; // something is wrong and/or unexpected

        struct string len_str, addr_str;
        string_init(&len_str);
        string_init(&addr_str);

        string_substr(&len_str, &dat_local, last_comma_idx + 1,
                      string_length(&dat_local) - 1);
        string_substr(&addr_str, &dat_local,
                      first_comma_idx + 1, last_comma_idx - 1);

        uint32_t length = string_read_hex32(&len_str, 0);
        addr32_t watch_addr = string_read_hex32(&addr_str, 0);;

        string_cleanup(&addr_str);
        string_cleanup(&len_str);

        int err_code = debug_add_w_watch(stub->dbg, watch_addr, length);

        if (err_code == 0)
            string_set(out, "OK");
        else
            err_str(out, err_code);

        goto cleanup;
    } else if (dat_c_str[1] == '3') {
        // read watchpoint

        if (string_find_first_of(&dat_local, ":") >= 0)
            goto cleanup; // we don't support conditions

        int first_comma_idx = string_find_first_of(&dat_local, ",");
        if ((first_comma_idx < 0) ||
            (first_comma_idx == (int)string_length(&dat_local) - 1)) {
            goto cleanup; // something is wrong and/or unexpected
        }

        int last_comma_idx = string_find_last_of(&dat_local, ",");
        if (last_comma_idx < 0)
            goto cleanup; // something is wrong and/or unexpected

        struct string len_str, addr_str;
        string_init(&len_str);
        string_init(&addr_str);

        string_substr(&len_str, &dat_local,
                      last_comma_idx + 1, string_length(&dat_local) - 1);
        string_substr(&addr_str, &dat_local,
                      first_comma_idx + 1, last_comma_idx - 1);

        uint32_t length = string_read_hex32(&len_str, 0);
        addr32_t watch_addr = string_read_hex32(&addr_str, 0);;

        string_cleanup(&addr_str);
        string_cleanup(&len_str);

        int err_code = debug_add_r_watch(stub->dbg, watch_addr, length);

        if (err_code == 0)
            string_set(out, "OK");
        else
            err_str(out, err_code);

        goto cleanup;
    } else {
        // unsupported
        goto cleanup;
    }

cleanup:
    string_cleanup(&dat_local);
}

static void handle_z_packet(struct gdb_stub *stub, struct string *out,
                            struct string const *dat) {
    char const *dat_c_str = string_get(dat);
    struct string dat_local;
    string_init(&dat_local);
    string_copy(&dat_local, dat);

    if (dat_c_str[1] == '1') {
        //hardware breakpoint

        if (string_find_first_of(&dat_local, ":") >= 0)
            goto cleanup; // we don't support conditions

        int first_comma_idx = string_find_first_of(&dat_local, ",");
        if ((first_comma_idx < 0) ||
            (first_comma_idx == (int)string_length(&dat_local) - 1)) {
            goto cleanup; // something is wrong and/or unexpected
        }

        int last_comma_idx = string_find_last_of(&dat_local, ",");
        if (last_comma_idx < 0)
            goto cleanup; // something is wrong and/or unexpected

        struct string tmp_str;
        string_init(&tmp_str);
        string_substr(&tmp_str, &dat_local,
                      first_comma_idx + 1, last_comma_idx - 1);
        string_cleanup(&dat_local);
        memcpy(&dat_local, &tmp_str, sizeof(dat_local));

        addr32_t break_addr = string_read_hex32(&dat_local, 0);

        int err_code = debug_remove_break(stub->dbg, break_addr);

        if (err_code == 0)
            string_set(out, "OK");
        else
            err_str(out, err_code);

        goto cleanup;
    } else if (dat_c_str[1] == '2') {
        // write watchpoint

        if (string_find_first_of(&dat_local, ":") >= 0)
            goto cleanup; // we don't support conditions

        int first_comma_idx = string_find_first_of(&dat_local, ",");
        if ((first_comma_idx < 0) ||
            (first_comma_idx == (int)string_length(&dat_local) - 1)) {
            goto cleanup; // something is wrong and/or unexpected
        }

        int last_comma_idx = string_find_last_of(&dat_local, ",");
        if (last_comma_idx < 0)
            goto cleanup; // something is wrong and/or unexpected

        struct string len_str, addr_str;
        string_init(&len_str);
        string_init(&addr_str);

        string_substr(&len_str, &dat_local, last_comma_idx + 1,
                      string_length(&dat_local) - 1);
        string_substr(&addr_str, &dat_local, first_comma_idx + 1, last_comma_idx - 1);

        uint32_t length = string_read_hex32(&len_str, 0);
        addr32_t watch_addr = string_read_hex32(&addr_str, 0);;

        string_cleanup(&addr_str);
        string_cleanup(&len_str);

        int err_code = debug_remove_w_watch(stub->dbg, watch_addr, length);

        if (err_code == 0)
            string_set(out, "OK");
        else
            err_str(out, err_code);

        goto cleanup;
    } else if (dat_c_str[1] == '3') {
        // read watchpoint

        if (string_find_first_of(&dat_local, ":") >= 0)
            goto cleanup; // we don't support conditions

        int first_comma_idx = string_find_first_of(&dat_local, ",");
        if ((first_comma_idx < 0) ||
            (first_comma_idx == (int)string_length(&dat_local) - 1)) {
            goto cleanup; // something is wrong and/or unexpected
        }

        int last_comma_idx = string_find_last_of(&dat_local, ",");
        if (last_comma_idx < 0)
            goto cleanup; // something is wrong and/or unexpected

        struct string len_str, addr_str;
        string_init(&len_str);
        string_init(&addr_str);

        string_substr(&len_str, &dat_local, last_comma_idx + 1,
                      string_length(&dat_local) - 1);
        string_substr(&addr_str, &dat_local,
                      first_comma_idx + 1, last_comma_idx - 1);

        uint32_t length = string_read_hex32(&len_str, 0);
        addr32_t watch_addr = string_read_hex32(&addr_str, 0);;

        string_cleanup(&addr_str);
        string_cleanup(&len_str);

        int err_code = debug_remove_r_watch(stub->dbg, watch_addr, length);

        if (err_code == 0)
            string_set(out, "OK");
        else
            err_str(out, err_code);

        goto cleanup;
    } else {
        // unsupported
        goto cleanup;
    }

cleanup:
    string_cleanup(&dat_local);
}

static void handle_packet(struct gdb_stub *stub,  struct string *pkt) {
    struct string dat;
    struct string response;
    struct string resp_pkt;

    string_init(&dat);
    string_init(&response);
    string_init(&resp_pkt);

    extract_packet(&dat, pkt);

    if (string_length(&dat)) {
        char first_ch = string_get(&dat)[0];
        if (first_ch == 'q') {
            handle_q_packet(stub, &response, &dat);
        } else if (first_ch == 'g') {
            handle_g_packet(stub, &response, &dat);
        } else if (first_ch == 'G') {
            handle_G_packet(stub, &response, &dat);
        } else if (first_ch == 'm') {
            handle_m_packet(stub, &response, &dat);
        } else if (first_ch == 'M') {
            handle_M_packet(stub, &response, &dat);
        } else if (first_ch == '?') {
            string_set(&response, "S05 create:");
        } else if (first_ch == 's') {
            handle_s_packet(stub, &response, &dat);
            goto cleanup;
        } else if (first_ch == 'c') {
            handle_c_packet(stub, &response, &dat);
            goto cleanup;
        } else if (first_ch == 'P') {
            handle_P_packet(stub, &response, &dat);
        } else if (first_ch == 'D') {
            handle_D_packet(stub, &response, &dat);
        } else if (first_ch == 'k') {
            handle_K_packet(stub, &response, &dat);
        } else if (first_ch == 'z') {
            handle_z_packet(stub, &response, &dat);
        } else if (first_ch == 'Z') {
            handle_Z_packet(stub, &response, &dat);
        }
    }

    craft_packet(&resp_pkt, &response);
    transmit_pkt(stub, &resp_pkt);

cleanup:
    string_cleanup(&resp_pkt);
    string_cleanup(&response);
    string_cleanup(&dat);
}

static void transmit_pkt(struct gdb_stub *stub, struct string const *pkt) {
#ifdef GDBSTUB_VERBOSE
    std::cout << ">>>> " << string_get(pkt) << std::endl;
#endif

    string_copy(&stub->unack_packet, pkt);
    transmit(stub, pkt);
}

static bool next_packet(struct gdb_stub *stub, struct string *pkt) {
    struct string pktbuf_tmp;
    struct string tmp_str;
    char ch;
    bool found_pkt = false;

    string_init(&pktbuf_tmp);
    string_copy(&pktbuf_tmp, &stub->input_packet);

    // wait around for the start character, ignore all other characters
    do {
        if (!string_length(&pktbuf_tmp))
            goto cleanup;

        ch = string_get(&pktbuf_tmp)[0];

        string_init(&tmp_str);
        string_substr(&tmp_str, &pktbuf_tmp, 1, string_length(&pktbuf_tmp) - 1);
        string_cleanup(&pktbuf_tmp);
        memcpy(&pktbuf_tmp, &tmp_str, sizeof(pktbuf_tmp));

        if (ch == std::char_traits<char>::eof())
            goto cleanup;
    } while(ch != '$');
    string_append_char(pkt, ch);

    // now, read until a # or end of buffer is found
    do {
        if (!string_length(&pktbuf_tmp))
            goto cleanup;

        ch = string_get(&pktbuf_tmp)[0];

        string_init(&tmp_str);
        string_substr(&tmp_str, &pktbuf_tmp, 1, string_length(&pktbuf_tmp) - 1);
        string_cleanup(&pktbuf_tmp);
        memcpy(&pktbuf_tmp, &tmp_str, sizeof(pktbuf_tmp));

        if (ch == std::char_traits<char>::eof())
            goto cleanup;

        string_append_char(pkt, ch);
    } while(ch != '#');
        
    // read the one-byte (two char) checksum
    if (!string_length(&pktbuf_tmp))
            goto cleanup;
    ch = string_get(&pktbuf_tmp)[0];
    if (ch == std::char_traits<char>::eof())
            goto cleanup;
    string_init(&tmp_str);
    string_substr(&tmp_str, &pktbuf_tmp, 1, string_length(&pktbuf_tmp) - 1);
    string_cleanup(&pktbuf_tmp);
    memcpy(&pktbuf_tmp, &tmp_str, sizeof(pktbuf_tmp));
    string_append_char(pkt, ch);

    if (!string_length(&pktbuf_tmp))
            goto cleanup;
    ch = string_get(&pktbuf_tmp)[0];
    if (ch == std::char_traits<char>::eof())
            goto cleanup;
    string_init(&tmp_str);
    string_substr(&tmp_str, &pktbuf_tmp, 1, string_length(&pktbuf_tmp) - 1);
    string_cleanup(&pktbuf_tmp);
    memcpy(&pktbuf_tmp, &tmp_str, sizeof(pktbuf_tmp));
    string_append_char(pkt, ch);

    string_set(&stub->input_packet, string_get(&pktbuf_tmp));

#ifdef GDBSTUB_VERBOSE
    std::cout << "<<<< " << string_get(pkt) << std::endl;
#endif

    found_pkt = true;

cleanup:
    string_cleanup(&pktbuf_tmp);
    return found_pkt;
}

// static void errror_handler(int error_tp, void *argptr) {
//     struct gdb_stub* stub = (struct gdb_stub*)argptr;

//     if (stub->should_expect_mem_access_error) {
//         stub->mem_access_error = true;
//     } else {
//         error_print();
//         exit(1);
//     }
// }

static void expect_mem_access_error(struct gdb_stub *stub, bool should) {
    stub->mem_access_error = false;
    stub->should_expect_mem_access_error = true;
}

static void
listener_cb(struct evconnlistener *listener,
            evutil_socket_t fd, struct sockaddr *saddr,
            int socklen, void *arg) {
    struct gdb_stub *stub = (struct gdb_stub*)arg;

    if (!stub->is_listening)
        return;

    stub->bev = bufferevent_socket_new(dc_event_base, fd,
                                       BEV_OPT_CLOSE_ON_FREE);
    if (!stub->bev)
        RAISE_ERROR(ERROR_FAILED_ALLOC);

    bufferevent_setcb(stub->bev, handle_read,
                      handle_write, handle_events, stub);
    bufferevent_enable(stub->bev, EV_WRITE);
    bufferevent_enable(stub->bev, EV_READ);

    stub->is_listening = false;
}

static void handle_events(struct bufferevent *bev, short events, void *arg) {
    exit(2);
}

static void handle_read(struct bufferevent *bev, void *arg) {
    struct evbuffer *read_buffer;
    struct gdb_stub *stub = (struct gdb_stub*)arg;

    if (!(read_buffer = evbuffer_new()))
        RAISE_ERROR(ERROR_FAILED_ALLOC);

    bufferevent_read_buffer(bev, read_buffer);
    size_t buflen = evbuffer_get_length(read_buffer);

    for (unsigned i = 0; i < buflen; i++) {
        uint8_t tmp;
        if (evbuffer_remove(read_buffer, &tmp, sizeof(tmp)) < 0)
            RAISE_ERROR(ERROR_FAILED_ALLOC);

        char c = tmp;

        if (string_length(&stub->input_packet)) {
            string_append_char(&stub->input_packet, c);

            struct string pkt;
            string_init(&pkt);
            bool pkt_valid = next_packet(stub, &pkt);
            if (string_length(&pkt) && pkt_valid) {
                string_set(&stub->input_packet, "");

                // TODO: verify the checksum

#ifdef GDBSTUB_VERBOSE
                std::cout << ">>>> +" << std::endl;
#endif
                struct string plus_symbol;
                string_init_txt(&plus_symbol, "+");
                transmit(stub, &plus_symbol);
                string_cleanup(&plus_symbol);
                handle_packet(stub, &pkt);
            }

            string_cleanup(&pkt);
        } else {
            if (c == '+') {
#ifdef GDBSTUB_VERBOSE
                std::cout << "<<<< +" << std::endl;
#endif
                if (!string_length(&stub->unack_packet))
                    std::cerr << "WARNING: received acknowledgement for unsent " <<
                        "packet" << std::endl;
                string_set(&stub->unack_packet, "");
            } else if (c == '-') {
#ifdef GDBSTUB_VERBOSE
                std::cout << "<<<< -" << std::endl;
#endif
                if (!string_length(&stub->unack_packet)) {
                    std::cerr << "WARNING: received negative acknowledgement for " <<
                        "unsent packet" << std::endl;
                } else {
#ifdef GDBSTUB_VERBOSE
                    std::cout << ">>>>" << string_get(&stub->unack_packet) << std::endl;
#endif
                    transmit(stub, &stub->unack_packet);
                }
            } else if (c == '$') {
                // new packet
                string_set(&stub->input_packet, "$");
            } else if (c == 3) {
                // user pressed ctrl+c (^C) on the gdb frontend
                std::cout << "GDBSTUB: user requested breakpoint (ctrl-C)" <<
                    std::endl;
                if (stub->dbg->cur_state == DEBUG_STATE_NORM) {
                    gdb_on_break(stub);
                    stub->dbg->cur_state = DEBUG_STATE_BREAK;
                }
            } else {
                std::cerr << "WARNING: ignoring unexpected character " << c <<
                    std::endl;
            }
        }
    }
}

/*
 * this function gets called when libevent is done writing
 * and is hungry for more data
 */
static void handle_write(struct bufferevent *bev, void *arg) {
    struct gdb_stub *stub = (struct gdb_stub*)arg;

    stub->is_writing = false;

    if (!evbuffer_get_length(stub->output_buffer))
        return;

    write_start(stub);
}
