/*******************************************************************************
 *
 *
 *    WashingtonDC Dreamcast Emulator
 *    Copyright (C) 2016 snickerbockers
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

#include <unistd.h>
#include <iostream>
#include <sstream>
#include <limits>

#include "BaseException.hpp"
#include "hw/sh4/Memory.hpp"
#include "hw/sh4/sh4.hpp"
#include "hw/sh4/Icache.hpp"
#include "hw/sh4/Ocache.hpp"
#include "tool/sh4asm/sh4asm.hpp"
#include "RandGenerator.hpp"

typedef RandGenerator<boost::uint32_t> RandGen32;
typedef int(*inst_test_func_t)(Sh4 *cpu, Memory *mem, RandGen32 *randgen32);

class Sh4InstTests {
public:
    /*
     * Put the cpu in a "clean" default state.
     */
    static void reset_cpu(Sh4 *cpu) {
        cpu->reg.pc = 0;
        cpu->reg.sr = Sh4::SR_MD_MASK;

#ifdef ENABLE_SH4_OCACHE
        cpu->op_cache->reset();
#endif
#ifdef ENABLE_SH4_ICACHE
        cpu->inst_cache->reset();
#endif
    }

    class Operand {
    public:
        Operand(Sh4 *cpu_ptr, Memory *mem_ptr) {
            this->cpu_ptr = cpu_ptr;
            this->mem_ptr = mem_ptr;
        }

        Sh4 *cpu() {
            return cpu_ptr;
        }

        Memory *mem() {
            return mem_ptr;
        }

    private:
        Sh4 *cpu_ptr;
        Memory *mem_ptr;
    };

    // An argument that consists of a single general-purpose register
    class GenRegArg : private Operand {
    public:
        GenRegArg(Sh4 *cpu, Memory *mem, unsigned reg_no, reg32_t reg_val) :
            Operand(cpu, mem){
            this->reg_no = reg_no;
            this->initial_val = this->final_val = reg_val;
        }

        GenRegArg(Sh4 *cpu, Memory *mem, unsigned reg_no,
                  reg32_t initial_val, reg32_t final_val) : Operand(cpu, mem) {
            this->reg_no = reg_no;
            this->initial_val = initial_val;
            this->final_val = final_val;
        }

        void initialize() {
            *(cpu()->gen_reg(reg_no)) = initial_val;
        }

        bool verify() {
            if (*(cpu()->gen_reg(reg_no)) == final_val)
                return true;
            std::cout << "ERROR ON REGISTER " << txt() << std::endl;
            std::cout << "EXPECTED VALUE: " << std::hex << final_val <<
                std::endl;
            std::cout << "ACTUAL VALUE: " << *(cpu()->gen_reg(reg_no)) <<
                std::endl;
            return false;
        }

        std::string txt() {
            std::stringstream ss;
            ss << "R" << reg_no;
            return ss.str();
        }

        reg32_t get_val() {
            return *cpu()->gen_reg(reg_no);
        }
    private:
        unsigned reg_no;
        reg32_t initial_val, final_val;
    };

    class SpecRegArg : private Operand {
    public:
        SpecRegArg(Sh4 *cpu, Memory *mem, char const *name,
                   reg32_t *ptr, reg32_t reg_val) : Operand(cpu, mem) {
            this->name = name;
            this->ptr = ptr;
            this->initial_val = this->final_val = reg_val;
        }

        SpecRegArg(Sh4 *cpu, Memory *mem, char const *name,
                   reg32_t *ptr, reg32_t initial_val, reg32_t final_val) :
            Operand(cpu, mem) {
            this->name = name;
            this->ptr = ptr;
            this->initial_val = initial_val;
            this->final_val = final_val;
        }

        void initialize() {
            *ptr = initial_val;
        }

        bool verify() {
            if (*ptr == final_val)
                return true;
            std::cout << "ERROR ON REGISTER " << txt() << std::endl;
            std::cout << "EXPECTED VALUE: " << std::hex << final_val <<
                std::endl;
            std::cout << "ACTUAL VALUE: " << *ptr << std::endl;
            return false;
        }

        std::string txt() {
            return std::string(name);
        }

        reg32_t get_val() {
            return *ptr;
        }
    private:
        char const *name;
        reg32_t *ptr;
        reg32_t initial_val, final_val;
    };

    class ImmedArg : private Operand {
    public:
        ImmedArg(Sh4 *cpu, Memory *mem, uint16_t val) : Operand(cpu, mem) {
            this->val = val;
        }

        void initialize() {
            // do nothing
        }

        bool verify() {
            return true;
        }

        std::string txt() {
            std::stringstream ss;
            ss << "#" << val;
            return ss.str();
        }

        uint16_t get_val() {
            return val;
        }
    private:
        uint16_t val;
    };

    class DispArg : private Operand {
    public:
        DispArg(Sh4 *cpu, Memory *mem, unsigned val) : Operand(cpu, mem) {
            this->val = val;
        }

        void initialize() {
            // do nothing
        }

        bool verify() {
            return true;
        }

        std::string txt() {
            std::stringstream ss;
            ss << val;
            return ss.str();
        }

        uint16_t get_val() {
            return val;
        }
    private:
        unsigned val;
    };

    template <class OpLeft, class OpRight, unsigned TargetSize>
    class BinIndArg : private Operand {
    public:
        /*
         * The cache_addr param is for situations where the value of an operand
         * might change but the old addr should be used for verifications (like
         * PC and @Rm+)
         */
        BinIndArg(Sh4 *cpu, Memory *mem, OpLeft lhs_op, OpRight rhs_op,
                  uint32_t val_init, uint32_t val_expect, unsigned lhs_scale = 1,
                  unsigned rhs_scale = 1, unsigned offset = 0,
                  bool cache_addr = false,
                  reg32_t lhs_addr_mask = ~0, reg32_t rhs_addr_mask = ~0) :
            Operand(cpu, mem), lhs(lhs_op), rhs(rhs_op) {
            this->lhs_scale = lhs_scale;
            this->rhs_scale = rhs_scale;
            this->offset = offset;
            this->val_init = val_init;
            this->val_expect = val_expect;
            this->cache_addr = cache_addr;
            this->lhs_addr_mask = lhs_addr_mask;
            this->rhs_addr_mask = rhs_addr_mask;
            this->cached_addr_valid = false;
        }

        void initialize() {
            lhs.initialize();
            rhs.initialize();
            cached_addr = get_addr();
            cpu()->write_mem(&val_init, cached_addr, TargetSize);
        }

        std::string txt() {
            std::stringstream ss;
            ss << "@(" << lhs.txt() << ", " << rhs.txt() << ")";
            return ss.str();
        }

        bool verify() {
            uint32_t val_actual = get_val();
            if (val_actual == val_expect)
                return true;

            std::cout << "Expected value was " << std::hex << val_expect <<
                std::endl << "Actual value is " << val_actual << std::endl <<
                "address is " << get_addr() << std::endl;
            return false;
        }

        addr32_t get_addr() {
            if (cache_addr && cached_addr_valid) {
                return cached_addr;
            }
            cached_addr_valid = true;
            cached_addr = addr32_t(lhs_addr_mask & lhs.get_val()) * lhs_scale +
                addr32_t(rhs_addr_mask & rhs.get_val()) * rhs_scale + offset;
            return cached_addr;
        }

        uint32_t get_val() {
            uint32_t val;
            addr32_t addr = get_addr();
            cpu()->read_mem(&val, addr, TargetSize);
            return val;
        }
    private:
        bool cache_addr, cached_addr_valid;
        addr32_t cached_addr;
        OpLeft lhs;
        OpRight rhs;
        uint32_t val_expect, val_init;
        unsigned lhs_scale, rhs_scale, offset;
        reg32_t lhs_addr_mask, rhs_addr_mask;
    };

    class DefValidationFunc {
    public:
        bool operator()(Sh4 *cpu, Memory *mem) const {
            return true;
        }
    };

    class ValidateFlagT {
    public:
        ValidateFlagT(bool expect) {
            this->expect = expect;
        }
        bool operator()(Sh4 *cpu, Memory *mem) const {
            bool tval = bool(!!(cpu->reg.sr & Sh4::SR_FLAG_T_MASK));
            bool res = tval == expect;

            if (!res)
                std::cout << "ERROR: T flag does not match (expected T=" <<
                    expect << ", got T=" << tval << ")" << std::endl;
            return res;
        }
    private:
        bool expect;
    };

    // common test infrastructure for binary unit_tests
    template <class SrcOpClass, class DstOpClass, class Validator>
    static int do_binary_reg_reg(Sh4 *cpu, Memory *mem,
                                 std::string opcode,
                                 SrcOpClass src_op, DstOpClass dst_op,
                                 Validator validate, addr32_t inst_loc = 0) {
        Sh4Prog test_prog;
        std::stringstream ss;
        std::string cmd;

        ss << opcode << " " << src_op.txt() << ", " << dst_op.txt() << "\n";
        cmd = ss.str();
        test_prog.assemble(cmd);
        const Sh4Prog::InstList& inst = test_prog.get_prog();
        mem->load_program(inst_loc, inst.begin(), inst.end());

        reset_cpu(cpu);

        src_op.initialize();
        dst_op.initialize();
        cpu->exec_inst();

        if (!src_op.verify() || !dst_op.verify() || !validate(cpu, mem)) {
            std::cout << "ERROR while running " << cmd << std::endl;
            reset_cpu(cpu); // in case SR (or other state register) changed
            return 1;
        }

        reset_cpu(cpu); // in case SR (or other state register) changed
        return 0;
    }

    // very basic test that does a whole lot of nothing
    static int nop_test(Sh4 *cpu, Memory *mem, RandGen32 *randgen32) {
        Sh4Prog test_prog;
        test_prog.assemble("NOP\n");
        const Sh4Prog::InstList& inst = test_prog.get_prog();
        mem->load_program(0, inst.begin(), inst.end());

        reset_cpu(cpu);

        cpu->exec_inst();

        return 0;
    }

    // ADD #imm, Rn
    // 0111nnnniiiiiiii
    static int add_immed_test(Sh4 *cpu, Memory *mem, RandGen32 *randgen32) {
        reg32_t initial_val = randgen32->pick_val(0);
        int failure = 0;
        for (int reg_no = 0; reg_no <= 15; reg_no++) {
            for (unsigned imm_val = 0; imm_val <= 0xff; imm_val++) {
                uint32_t dst_val = randgen32->pick_val(0);
                ImmedArg src_op(cpu, mem, imm_val);
                GenRegArg dst_op(cpu, mem, reg_no, dst_val,
                                 int32_t(dst_val) + int32_t(int8_t(imm_val)));
                failure = failure ||
                    do_binary_reg_reg(cpu, mem, "ADD", src_op, dst_op,
                                      DefValidationFunc());
            }
        }
        return failure;
    }

    // ADD Rm, Rn
    // 0111nnnnmmmm1100
    static int add_gen_gen_test(Sh4 *cpu, Memory *mem, RandGen32 *randgen32) {
        int failure = 0;
        for (int reg1_no = 0; reg1_no <= 15; reg1_no++) {
            for (int reg2_no = 0; reg2_no <= 15; reg2_no++) {
                reg32_t src_val = randgen32->pick_val(0);
                reg32_t dst_val = (reg1_no == reg2_no ? src_val :
                                   randgen32->pick_val(0));
                GenRegArg src(cpu, mem, reg1_no, src_val,
                              src_val == dst_val ? src_val + dst_val : src_val);
                GenRegArg dst(cpu, mem, reg2_no, dst_val, src_val + dst_val);
                failure = failure || do_binary_reg_reg(cpu, mem, "ADD",
                                                       src, dst,
                                                       DefValidationFunc());
            }
        }
        return failure;
    }

    // ADDC Rm, Rn
    // 0011nnnnmmmm1110
    static int do_addc_gen_gen_test(Sh4 *cpu, Memory *mem,
                                    reg32_t src1, reg32_t src2) {
        int failure = 0;

        for (int reg1_no = 0; reg1_no <= 15; reg1_no++) {
            for (int reg2_no = 0; reg2_no <= 15; reg2_no++) {
                bool expected_t;
                uint32_t initial_val1 = src1;
                uint32_t initial_val2;

                if (reg1_no == reg2_no)
                    initial_val2 = initial_val1;
                else
                    initial_val2 = src2;

                int32_t sum = initial_val1 + initial_val2;
                GenRegArg reg1(cpu, mem, reg1_no, initial_val1,
                               reg1_no == reg2_no ? sum : initial_val1);
                GenRegArg reg2(cpu, mem, reg2_no, initial_val2, sum);

                bool carry;
                uint64_t expected_val64 = uint64_t(initial_val1) +
                    uint64_t(initial_val2);
                if (expected_val64 == uint64_t(initial_val1 + initial_val2)) {
                    carry = false;
                } else {
                    carry = true;
                }

                failure = failure ||
                    do_binary_reg_reg(cpu, mem, "ADDC", reg1, reg2,
                                      ValidateFlagT(carry));
            }
        }
        return failure;
    }

    // ADDC Rm, Rn
    // 0011nnnnmmmm1110
    static int addc_gen_gen_test(Sh4 *cpu, Memory *mem, RandGen32 *randgen32) {
        int failed = 0;

        // run the test with a couple random values
        failed = failed || do_addc_gen_gen_test(cpu, mem,
                                                randgen32->pick_val(0),
                                                randgen32->pick_val(0));

        // make sure we get at least one value in that should not cause a carry
        failed = failed || do_addc_gen_gen_test(cpu, mem, 0, 0);

        // make sure we get at least one value in that should cause a carry
        failed = failed ||
            do_addc_gen_gen_test(cpu, mem, std::numeric_limits<reg32_t>::max(),
                                 std::numeric_limits<reg32_t>::max());

        // test a value that should *almost* cause a carry
        failed = failed ||
            do_addc_gen_gen_test(cpu, mem, 1,
                                 std::numeric_limits<reg32_t>::max() - 1);

        // test a value pair that should barely cause a carry
        failed = failed ||
            do_addc_gen_gen_test(cpu, mem,
                                 std::numeric_limits<reg32_t>::max() - 1, 2);

        return failed;
    }

    // ADDV Rm, Rn
    // 0011nnnnmmmm1111
    static int do_addv_gen_gen_test(Sh4 *cpu, Memory *mem,
                                    reg32_t src1, reg32_t src2) {
        int failure = 0;

        for (int reg1_no = 0; reg1_no <= 15; reg1_no++) {
            for (int reg2_no = 0; reg2_no <= 15; reg2_no++) {
                bool expected_t;
                int32_t initial_val1 = src1;
                int32_t initial_val2;

                if (reg1_no == reg2_no)
                    initial_val2 = initial_val1;
                else
                    initial_val2 = src2;

                int32_t sum = initial_val1 + initial_val2;
                GenRegArg reg1(cpu, mem, reg1_no, initial_val1,
                               reg1_no == reg2_no ? sum : initial_val1);
                GenRegArg reg2(cpu, mem, reg2_no, initial_val2, sum);

                uint64_t expected_val64 = uint64_t(uint32_t(initial_val1)) +
                    uint64_t(uint32_t(initial_val2));

                bool overflow;

                if (initial_val1 >= 0 && initial_val2 >= 0) {
                    overflow = ((std::numeric_limits<int32_t>::max() -
                                 initial_val1) < initial_val2);
                } else {
                    overflow = ((std::numeric_limits<int32_t>::min() -
                                 initial_val2) > initial_val1);
                }

                failure = failure ||
                    do_binary_reg_reg(cpu, mem, "ADDV", reg1, reg2,
                                      ValidateFlagT(overflow));
            }
        }
        return failure;
    }

    // ADDV Rm, Rn
    // 0011nnnnmmmm1111
    static int addv_gen_gen_test(Sh4 *cpu, Memory *mem, RandGen32 *randgen32) {
        int failed = 0;
        randgen32->reset();

        // this should not overflow
        failed = failed || do_addv_gen_gen_test(cpu, mem, 0, 0);

        // random values for good measure
        failed = failed || do_addv_gen_gen_test(cpu, mem,
                                                randgen32->pick_val(0),
                                                randgen32->pick_val(0));

        // *almost* overflow positive to negative
        failed = failed ||
            do_addv_gen_gen_test(cpu, mem, 1,
                                 std::numeric_limits<int32_t>::max() - 1);

        // slight overflow positive to negative
        failed = failed ||
            do_addv_gen_gen_test(cpu, mem, 2,
                                 std::numeric_limits<int32_t>::max() - 1);

        // massive overflow positive to negative
        failed = failed ||
            do_addv_gen_gen_test(cpu, mem,
                                 std::numeric_limits<int32_t>::max(),
                                 std::numeric_limits<int32_t>::max());

        // *almost* overflow negative to positive
        failed = failed ||
            do_addv_gen_gen_test(cpu, mem,
                                 std::numeric_limits<int32_t>::min() + 1, 1);

        // slight overflow negative to positive
        failed = failed ||
            do_addv_gen_gen_test(cpu, mem,
                                 std::numeric_limits<int32_t>::min() + 1, 2);

        // massive overflow negative to positive
        failed = failed ||
            do_addv_gen_gen_test(cpu, mem,
                                 std::numeric_limits<int32_t>::min(),
                                 std::numeric_limits<int32_t>::min());

        return failed;
    }

    // SUB Rm, Rn
    // 0011nnnnmmmm1000
    static int sub_gen_gen_test(Sh4 *cpu, Memory *mem, RandGen32 *randgen32) {
        int failure = 0;
        for (int reg1_no = 0; reg1_no <= 15; reg1_no++) {
            for (int reg2_no = 0; reg2_no <= 15; reg2_no++) {
                reg32_t src_val = randgen32->pick_val(0);
                reg32_t dst_val = (reg1_no == reg2_no ? src_val :
                                   randgen32->pick_val(0));
                GenRegArg src(cpu, mem, reg1_no, src_val,
                              src_val == dst_val ? dst_val - src_val : src_val);
                GenRegArg dst(cpu, mem, reg2_no, dst_val, dst_val - src_val);
                failure = failure || do_binary_reg_reg(cpu, mem, "SUB",
                                                       src, dst,
                                                       DefValidationFunc());
            }
        }
        return failure;
    }

    // SUBC Rm, Rn
    // 0011nnnnmmmm1010
    static int do_subc_gen_gen_test(Sh4 *cpu, Memory *mem,
                                    reg32_t src1, reg32_t src2) {
        int failure = 0;

        for (int reg1_no = 0; reg1_no <= 15; reg1_no++) {
            for (int reg2_no = 0; reg2_no <= 15; reg2_no++) {
                bool expected_t;
                uint32_t initial_val1 = src1;
                uint32_t initial_val2;

                if (reg1_no == reg2_no)
                    initial_val2 = initial_val1;
                else
                    initial_val2 = src2;

                int32_t diff = initial_val2 - initial_val1;
                GenRegArg reg1(cpu, mem, reg1_no, initial_val1,
                               reg1_no == reg2_no ? diff : initial_val1);
                GenRegArg reg2(cpu, mem, reg2_no, initial_val2, diff);

                bool carry = (initial_val1 > initial_val2);

                failure = failure ||
                    do_binary_reg_reg(cpu, mem, "SUBC", reg1, reg2,
                                      ValidateFlagT(carry));
            }
        }
        return failure;
    }

    static int subc_gen_gen_test(Sh4 *cpu, Memory *mem, RandGen32 *randgen32) {
        int failed = 0;

        // run the test with a couple random values
        failed = failed || do_subc_gen_gen_test(cpu, mem,
                                                randgen32->pick_val(0),
                                                randgen32->pick_val(0));

        // make sure we get at least one value in that should not cause a carry
        failed = failed || do_subc_gen_gen_test(cpu, mem, 0, 0);

        // make sure we get at least one value in that should cause a carry
        failed = failed ||
            do_subc_gen_gen_test(cpu, mem,
                                 std::numeric_limits<reg32_t>::max(), 0);

        // test a value that should *almost* cause a carry
        failed = failed ||
            do_subc_gen_gen_test(cpu, mem, std::numeric_limits<reg32_t>::max(),
                                 std::numeric_limits<reg32_t>::max());

        // test a value pair that should barely cause a carry
        failed = failed ||
            do_subc_gen_gen_test(cpu, mem, 1, 0);

        return failed;
    }

    static int do_subv_gen_gen_test(Sh4 *cpu, Memory *mem,
                                    reg32_t src1, reg32_t src2) {
        int failure = 0;

        for (int reg1_no = 0; reg1_no <= 15; reg1_no++) {
            for (int reg2_no = 0; reg2_no <= 15; reg2_no++) {
                bool expected_t;
                int32_t initial_val1 = src1;
                int32_t initial_val2;

                if (reg1_no == reg2_no)
                    initial_val2 = initial_val1;
                else
                    initial_val2 = src2;

                int32_t diff = initial_val2 - initial_val1;
                GenRegArg reg1(cpu, mem, reg1_no, initial_val1,
                               reg1_no == reg2_no ? diff : initial_val1);
                GenRegArg reg2(cpu, mem, reg2_no, initial_val2, diff);

                bool overflow = false;
                if (initial_val2 >= 0 && initial_val1 < 0) {
                    overflow = diff < 0;
                } else if (initial_val2 < 0 && initial_val1 >= 0) {
                    overflow = diff > 0;
                }

                failure = failure ||
                    do_binary_reg_reg(cpu, mem, "SUBV", reg1, reg2,
                                      ValidateFlagT(overflow));
            }
        }
        return failure;
    }

    static int subv_gen_gen_test(Sh4 *cpu, Memory *mem, RandGen32 *randgen32) {
        int failed = 0;

        // do one at random...
        failed = failed ||
            do_subv_gen_gen_test(cpu, mem,
                                 randgen32->pick_val(0), randgen32->pick_val(0));

        // now do one that's trivial
        failed = failed || do_subv_gen_gen_test(cpu, mem, 0, 0);

        // now do one that *almost* causes a negative overflow
        failed = failed ||
            do_subv_gen_gen_test(cpu, mem,
                                 -(std::numeric_limits<int32_t>::min() + 1), 0);

        // now do one that barely causes a negative overflow
        failed = failed ||
            do_subv_gen_gen_test(cpu, mem,
                                 -(std::numeric_limits<int32_t>::min() + 1), -1);

        // now do a massive negative overflow
        failed = failed ||
            do_subv_gen_gen_test(cpu, mem,
                                 -(std::numeric_limits<int32_t>::min() + 1),
                                 std::numeric_limits<int32_t>::min());

        // now do one that *almost* causes a positive overflow
        failed = failed ||
            do_subv_gen_gen_test(cpu, mem,
                                 -std::numeric_limits<int32_t>::max(), 0);

        // now do one that barely causes a positive overflow
        failed = failed ||
            do_subv_gen_gen_test(cpu, mem,
                                 -std::numeric_limits<int32_t>::max(), 1);

        // now do a massive positive overflow
        failed = failed ||
            do_subv_gen_gen_test(cpu, mem,
                                 -std::numeric_limits<int32_t>::max(),
                                 std::numeric_limits<int32_t>::max());

        return failed;
    }

    // MOVT Rn
    // 0000nnnn00101001
    static int movt_unary_gen_test(Sh4 *cpu, Memory *mem,
                                   RandGen32 *randgen32) {
        for (int reg_no = 0; reg_no < 16; reg_no++) {
            for (int t_val = 0; t_val < 2; t_val++) {
                Sh4Prog test_prog;
                std::stringstream ss;

                ss << "MOVT R" << reg_no << "\n";
                test_prog.assemble(ss.str());
                const Sh4Prog::InstList& inst = test_prog.get_prog();
                mem->load_program(0, inst.begin(), inst.end());

                reset_cpu(cpu);

                cpu->reg.sr &= ~Sh4::SR_FLAG_T_MASK;
                if (t_val)
                    cpu->reg.sr |= Sh4::SR_FLAG_T_MASK;

                 cpu->exec_inst();

                if (*cpu->gen_reg(reg_no) != t_val)
                    return 1;
            }
        }

        return 0;
    }

    // MOV #imm, Rn
    // 1110nnnniiiiiiii
    static int mov_binary_imm_gen_test(Sh4 *cpu, Memory *mem,
                                       RandGen32 *randgen32) {
        int failure = 0;
        for (int reg_no = 0; reg_no < 16; reg_no++) {
            reg32_t initial_val = randgen32->pick_val(0);
            for (unsigned imm_val = 0; imm_val <= 0xff; imm_val++) {
                uint32_t dst_val = randgen32->pick_val(0);
                ImmedArg src_op(cpu, mem, imm_val);
                GenRegArg dst_op(cpu, mem, reg_no, dst_val,
                                 int32_t(int8_t(imm_val)));
                failure = failure ||
                    do_binary_reg_reg(cpu, mem, "MOV", src_op, dst_op,
                                      DefValidationFunc());
            }
        }

        return failure;
    }

    // MOV.W @(disp, PC), Rn
    // 1001nnnndddddddd
    static int movw_binary_binind_disp_pc_gen(Sh4 *cpu, Memory *mem,
                                              RandGen32 *randgen32) {
        int failed = 0;

        for (int disp = 0; disp < 256; disp++) {
            for (int reg_no = 0; reg_no < 16; reg_no++) {
                addr32_t pc_max = mem->get_size() - 1 - 4 - disp * 2;
                addr32_t pc_val =
                    randgen32->pick_range(0, pc_max) & ~1;
                uint16_t val = randgen32->pick_val(0) & 0xffff;
                SpecRegArg pc_arg(cpu, mem, "PC", &cpu->reg.pc, pc_val);
                DispArg immed_arg(cpu, mem, disp);
                BinIndArg<DispArg, SpecRegArg, 2> src_op(
                    cpu, mem, immed_arg, pc_arg, val, val, 2, 1, 4, true);
                GenRegArg dst_op(cpu, mem, reg_no, randgen32->pick_val(0),
                                 int32_t(int16_t(val)));
                failed = failed ||
                    do_binary_reg_reg(cpu, mem, "MOV.W", src_op, dst_op,
                                      DefValidationFunc(), pc_val);
            }
        }
        return failed;
    }

    // MOV.L @(disp, PC), Rn
    // 1101nnnndddddddd
    static int movl_binary_binind_disp_pc_gen(Sh4 *cpu, Memory *mem,
                                              RandGen32 *randgen32) {
        int failed = 0;

        for (int disp = 0; disp < 256; disp++) {
            for (int reg_no = 0; reg_no < 16; reg_no++) {
                addr32_t pc_max = (mem->get_size() - 1 - 4 - disp * 4) | 3;
                addr32_t pc_val =
                    randgen32->pick_range(0, pc_max) & ~1;
                uint32_t val = randgen32->pick_val(0);
                SpecRegArg pc_arg(cpu, mem, "PC", &cpu->reg.pc, pc_val);
                DispArg immed_arg(cpu, mem, disp);
                BinIndArg<DispArg, SpecRegArg, 4> src_op(
                    cpu, mem, immed_arg, pc_arg, val, val, 4, 1, 4,
                    true, ~0, ~3);
                GenRegArg dst_op(cpu, mem, reg_no, randgen32->pick_val(0), val);
                failed = failed ||
                    do_binary_reg_reg(cpu, mem, "MOV.L", src_op, dst_op,
                                      DefValidationFunc(), pc_val);
            }
        }
        return failed;
    }

    // MOV Rm, Rn
    // 0110nnnnmmmm0011
    static int mov_binary_gen_gen(Sh4 *cpu, Memory *mem, RandGen32 *randgen32) {
        int failed = 0;

        for (unsigned reg_src = 0; reg_src < 16; reg_src++) {
            for (unsigned reg_dst = 0; reg_dst < 16; reg_dst++) {
                reg32_t src_val = randgen32->pick_val(0);
                reg32_t dst_val = (reg_src == reg_dst ?
                                   src_val : randgen32->pick_val(0));
                failed = failed ||
                    do_binary_reg_reg(cpu, mem, "MOV",
                                      GenRegArg(cpu, mem, reg_src, src_val),
                                      GenRegArg(cpu, mem, reg_dst,
                                                dst_val, src_val),
                                      DefValidationFunc());
            }
        }
        return failed;
    }

    // MOV.B Rm, @Rn
    // 0010nnnnmmmm0000
    static int do_movb_binary_gen_indgen(Sh4 *cpu, Memory *mem,
                                         unsigned addr, uint8_t val,
                                         unsigned reg_src, unsigned reg_dst) {
        Sh4Prog test_prog;
        std::stringstream ss;
        std::string cmd;
        uint8_t mem_val;

        if (reg_src == reg_dst)
            val = addr;

        ss << "MOV.B R" << reg_src << ", @R" << reg_dst << "\n";
        cmd = ss.str();
        test_prog.assemble(cmd);
        const Sh4Prog::InstList& inst = test_prog.get_prog();
        mem->load_program(0, inst.begin(), inst.end());

        reset_cpu(cpu);
        *cpu->gen_reg(reg_src) = val;
        *cpu->gen_reg(reg_dst) = addr;
        cpu->exec_inst();

        cpu->read_mem(&mem_val, addr, sizeof(mem_val));

        if (mem_val != val) {
            std::cout << "While running: " << cmd << std::endl;
            std::cout << "val is " << std::hex << (unsigned)val << std::endl;
            std::cout << "addr is " << std::hex << addr << std::endl;
            std::cout << "actual val is " << std::hex <<
                (unsigned)mem_val << std::endl;
            return 1;
        }

        return 0;
    }

    static int movb_binary_gen_indgen(Sh4 *cpu, Memory *mem,
                                      RandGen32 *randgen32) {
        int failed = 0;

        for (unsigned reg_src = 0; reg_src < 16; reg_src++) {
            for (unsigned reg_dst = 0; reg_dst < 16; reg_dst++) {
                addr32_t addr = randgen32->pick_val(0) % mem->get_size();
                failed = failed ||
                    do_movb_binary_gen_indgen(cpu, mem, addr,
                                              randgen32->pick_val(0) % 0xff,
                                              reg_src, reg_dst);
            }
        }

        return failed;
    }

    // MOV.W Rm, @Rn
    // 0010nnnnmmmm0001
    static int do_movw_binary_gen_indgen(Sh4 *cpu, Memory *mem,
                                         unsigned addr, uint16_t val,
                                         unsigned reg_src, unsigned reg_dst) {
        Sh4Prog test_prog;
        std::stringstream ss;
        std::string cmd;
        uint16_t mem_val;

        if (reg_src == reg_dst)
            val = addr;

        ss << "MOV.W R" << reg_src << ", @R" << reg_dst << "\n";
        cmd = ss.str();
        test_prog.assemble(cmd);
        const Sh4Prog::InstList& inst = test_prog.get_prog();
        mem->load_program(0, inst.begin(), inst.end());

        reset_cpu(cpu);
        *cpu->gen_reg(reg_src) = val;
        *cpu->gen_reg(reg_dst) = addr;
        cpu->exec_inst();

        cpu->read_mem(&mem_val, addr, sizeof(mem_val));

        if (mem_val != val) {
            std::cout << "While running: " << cmd << std::endl;
            std::cout << "val is " << std::hex << (unsigned)val << std::endl;
            std::cout << "addr is " << std::hex << addr << std::endl;
            std::cout << "actual val is " << std::hex <<
                (unsigned)mem_val << std::endl;
            return 1;
        }

        return 0;
    }

    static int movw_binary_gen_indgen(Sh4 *cpu, Memory *mem,
                                      RandGen32 *randgen32) {
        int failed = 0;

        for (unsigned reg_src = 0; reg_src < 16; reg_src++) {
            for (unsigned reg_dst = 0; reg_dst < 16; reg_dst++) {
                addr32_t addr = randgen32->pick_val(0) % (mem->get_size() - 1);
                failed = failed ||
                    do_movb_binary_gen_indgen(cpu, mem, addr,
                                              randgen32->pick_val(0) % 0xffff,
                                              reg_src, reg_dst);
            }
        }

        return failed;
    }

    // MOV.L Rm, @Rn
    // 0010nnnnmmmm0010
    static int do_movl_binary_gen_indgen(Sh4 *cpu, Memory *mem,
                                         unsigned addr, uint32_t val,
                                         unsigned reg_src, unsigned reg_dst) {
        Sh4Prog test_prog;
        std::stringstream ss;
        std::string cmd;
        uint8_t mem_val;

        if (reg_src == reg_dst)
            val = addr;

        ss << "MOV.L R" << reg_src << ", @R" << reg_dst << "\n";
        cmd = ss.str();
        test_prog.assemble(cmd);
        const Sh4Prog::InstList& inst = test_prog.get_prog();
        mem->load_program(0, inst.begin(), inst.end());

        reset_cpu(cpu);
        *cpu->gen_reg(reg_src) = val;
        *cpu->gen_reg(reg_dst) = addr;
        cpu->exec_inst();

        cpu->read_mem(&mem_val, addr, sizeof(mem_val));

        if (mem_val != val) {
            std::cout << "While running: " << cmd << std::endl;
            std::cout << "val is " << std::hex << (unsigned)val << std::endl;
            std::cout << "addr is " << std::hex << addr << std::endl;
            std::cout << "actual val is " << std::hex <<
                (unsigned)mem_val << std::endl;
            return 1;
        }

        return 0;
    }

    static int movl_binary_gen_indgen(Sh4 *cpu, Memory *mem,
                                      RandGen32 *randgen32) {
        int failed = 0;

        for (unsigned reg_src = 0; reg_src < 16; reg_src++) {
            for (unsigned reg_dst = 0; reg_dst < 16; reg_dst++) {
                addr32_t addr = randgen32->pick_val(0) % (mem->get_size() - 3);
                failed = failed ||
                    do_movb_binary_gen_indgen(cpu, mem, addr,
                                              randgen32->pick_val(0),
                                              reg_src, reg_dst);
            }
        }

        return failed;
    }

    // MOV.B @Rm, Rn
    // 0110nnnnmmmm0000
    static int do_movb_binary_indgen_gen(Sh4 *cpu, Memory *mem,
                                         unsigned addr, int8_t val,
                                         unsigned reg_src, unsigned reg_dst) {
        Sh4Prog test_prog;
        std::stringstream ss;
        std::string cmd;

        if (reg_src == reg_dst)
            val = addr;

        ss << "MOV.B @R" << reg_src << ", R" << reg_dst << "\n";
        cmd = ss.str();
        test_prog.assemble(cmd);
        const Sh4Prog::InstList& inst = test_prog.get_prog();
        mem->load_program(0, inst.begin(), inst.end());

        reset_cpu(cpu);
        *cpu->gen_reg(reg_src) = addr;
        cpu->write_mem(&val, addr, sizeof(val));
        cpu->exec_inst();

        if (*cpu->gen_reg(reg_dst) != int32_t(val)) {
            std::cout << "While running: " << cmd << std::endl;
            std::cout << "val is " << std::hex << (unsigned)val << std::endl;
            std::cout << "addr is " << std::hex << addr << std::endl;
            std::cout << "actual val is " << std::hex <<
                *cpu->gen_reg(reg_dst) << std::endl;
            return 1;
        }

        return 0;
    }

    static int movb_binary_indgen_gen(Sh4 *cpu, Memory *mem,
                                      RandGen32 *randgen32) {
        int failed = 0;

        for (unsigned reg_src = 0; reg_src < 16; reg_src++) {
            for (unsigned reg_dst = 0; reg_dst < 16; reg_dst++) {
                addr32_t addr = randgen32->pick_val(0) % mem->get_size();
                failed = failed ||
                    do_movb_binary_indgen_gen(cpu, mem, addr,
                                              randgen32->pick_val(0) % 0xff,
                                              reg_src, reg_dst);
            }
        }

        return failed;
    }

    // MOV.W @Rm, Rn
    // 0110nnnnmmmm0001
    static int do_movw_binary_indgen_gen(Sh4 *cpu, Memory *mem,
                                         unsigned addr, int16_t val,
                                         unsigned reg_src, unsigned reg_dst) {
        Sh4Prog test_prog;
        std::stringstream ss;
        std::string cmd;

        if (reg_src == reg_dst)
            val = addr;

        ss << "MOV.W @R" << reg_src << ", R" << reg_dst << "\n";
        cmd = ss.str();
        test_prog.assemble(cmd);
        const Sh4Prog::InstList& inst = test_prog.get_prog();
        mem->load_program(0, inst.begin(), inst.end());

        reset_cpu(cpu);
        *cpu->gen_reg(reg_src) = addr;
        cpu->write_mem(&val, addr, sizeof(val));
        cpu->exec_inst();

        if (*cpu->gen_reg(reg_dst) != int32_t(val)) {
            std::cout << "While running: " << cmd << std::endl;
            std::cout << "val is " << std::hex << (unsigned)val << std::endl;
            std::cout << "addr is " << std::hex << addr << std::endl;
            std::cout << "actual val is " << std::hex <<
                *cpu->gen_reg(reg_dst) << std::endl;
            return 1;
        }

        return 0;
    }

    static int movw_binary_indgen_gen(Sh4 *cpu, Memory *mem,
                                      RandGen32 *randgen32) {
        int failed = 0;

        for (unsigned reg_src = 0; reg_src < 16; reg_src++) {
            for (unsigned reg_dst = 0; reg_dst < 16; reg_dst++) {
                addr32_t addr = randgen32->pick_val(0) % (mem->get_size() - 1);
                failed = failed ||
                    do_movw_binary_indgen_gen(cpu, mem, addr,
                                              randgen32->pick_val(0) % 0xff,
                                              reg_src, reg_dst);
            }
        }

        return failed;
    }

    // MOV.L @Rm, Rn
    // 0110nnnnmmmm0010
    static int do_movl_binary_indgen_gen(Sh4 *cpu, Memory *mem,
                                         unsigned addr, int32_t val,
                                         unsigned reg_src, unsigned reg_dst) {
        Sh4Prog test_prog;
        std::stringstream ss;
        std::string cmd;

        if (reg_src == reg_dst)
            val = addr;

        ss << "MOV.L @R" << reg_src << ", R" << reg_dst << "\n";
        cmd = ss.str();
        test_prog.assemble(cmd);
        const Sh4Prog::InstList& inst = test_prog.get_prog();
        mem->load_program(0, inst.begin(), inst.end());

        reset_cpu(cpu);
        *cpu->gen_reg(reg_src) = addr;
        cpu->write_mem(&val, addr, sizeof(val));
        cpu->exec_inst();

        if (*cpu->gen_reg(reg_dst) != val) {
            std::cout << "While running: " << cmd << std::endl;
            std::cout << "val is " << std::hex << (unsigned)val << std::endl;
            std::cout << "addr is " << std::hex << addr << std::endl;
            std::cout << "actual val is " << std::hex <<
                *cpu->gen_reg(reg_dst) << std::endl;
            return 1;
        }

        return 0;
    }

    static int movl_binary_indgen_gen(Sh4 *cpu, Memory *mem,
                                      RandGen32 *randgen32) {
        int failed = 0;

        for (unsigned reg_src = 0; reg_src < 16; reg_src++) {
            for (unsigned reg_dst = 0; reg_dst < 16; reg_dst++) {
                addr32_t addr = randgen32->pick_val(0) % (mem->get_size() - 3);
                failed = failed ||
                    do_movw_binary_indgen_gen(cpu, mem, addr,
                                              randgen32->pick_val(0) % 0xff,
                                              reg_src, reg_dst);
            }
        }

        return failed;
    }

    // MOV.B Rm, @-Rn
    // 0010nnnnmmmm0100
    static int do_movb_binary_gen_inddecgen(Sh4 *cpu, Memory *mem,
                                            unsigned addr, uint8_t val,
                                            unsigned reg_src,
                                            unsigned reg_dst) {
        Sh4Prog test_prog;
        std::stringstream ss;
        std::string cmd;
        uint8_t mem_val;

        // increment addr 'cause the opcode is going to decrement it
        addr++;

        if (reg_src == reg_dst)
            val = addr;

        ss << "MOV.B R" << reg_src << ", @-R" << reg_dst << "\n";
        cmd = ss.str();
        test_prog.assemble(cmd);
        const Sh4Prog::InstList& inst = test_prog.get_prog();
        mem->load_program(0, inst.begin(), inst.end());

        reset_cpu(cpu);
        *cpu->gen_reg(reg_src) = val;
        *cpu->gen_reg(reg_dst) = addr;
        cpu->exec_inst();

        cpu->read_mem(&mem_val, addr-1, sizeof(mem_val));

        if (reg_src == reg_dst) {
            // special case - val will be decremented because the source and
            // destination are the same register
            val -= 1;
        }

        if (mem_val != val) {
            std::cout << "While running: " << cmd << std::endl;
            std::cout << "val is " << std::hex << (unsigned)val << std::endl;
            std::cout << "addr is " << std::hex << addr << std::endl;
            std::cout << "actual val is " << std::hex <<
                (unsigned)mem_val << std::endl;
            return 1;
        }

        if (*cpu->gen_reg(reg_dst) != addr - 1) {
            std::cout << "While running: " << cmd << std::endl;
            std::cout << "Expected the destination to be decremented "
                "(it was not)" << std::endl;
            std::cout << "val is " << std::hex << (unsigned)val << std::endl;
            std::cout << "addr is " << std::hex << addr << std::endl;
            std::cout << "actual val is " << std::hex <<
                (unsigned)mem_val << std::endl;
            return 1;
        }

        return 0;
    }

    static int movb_binary_gen_inddecgen(Sh4 *cpu, Memory *mem,
                                         RandGen32 *randgen32) {
        int failed = 0;

        for (int reg_src = 0; reg_src < 16; reg_src++) {
            for (int reg_dst = 0; reg_dst < 16; reg_dst++) {
                addr32_t addr = randgen32->pick_range(1, mem->get_size() - 2);
                failed = failed ||
                    do_movb_binary_gen_inddecgen(cpu, mem, addr,
                                                 randgen32->pick_val(0),
                                                 reg_src, reg_dst);
            }
        }

        return failed;
    }

    // MOV.W Rm, @-Rn
    // 0010nnnnmmmm0101
    static int do_movw_binary_gen_inddecgen(Sh4 *cpu, Memory *mem,
                                            unsigned addr, uint16_t val,
                                            unsigned reg_src,
                                            unsigned reg_dst) {
        Sh4Prog test_prog;
        std::stringstream ss;
        std::string cmd;
        uint16_t mem_val;

        // increment addr 'cause the opcode is going to decrement it
        addr += 2;

        if (reg_src == reg_dst)
            val = addr;

        ss << "MOV.W R" << reg_src << ", @-R" << reg_dst << "\n";
        cmd = ss.str();
        test_prog.assemble(cmd);
        const Sh4Prog::InstList& inst = test_prog.get_prog();
        mem->load_program(0, inst.begin(), inst.end());

        reset_cpu(cpu);
        *cpu->gen_reg(reg_src) = val;
        *cpu->gen_reg(reg_dst) = addr;
        cpu->exec_inst();

        cpu->read_mem(&mem_val, addr-2, sizeof(mem_val));

        if (reg_src == reg_dst) {
            // special case - val will be decremented because the source and
            // destination are the same register
            val -= 2;
        }

        if (mem_val != val) {
            std::cout << "While running: " << cmd << std::endl;
            std::cout << "val is " << std::hex << (unsigned)val << std::endl;
            std::cout << "addr is " << std::hex << addr << std::endl;
            std::cout << "actual val is " << std::hex <<
                (unsigned)mem_val << std::endl;
            return 1;
        }

        if (*cpu->gen_reg(reg_dst) != addr - 2) {
            std::cout << "While running: " << cmd << std::endl;
            std::cout << "Expected the destination to be decremented "
                "(it was not)" << std::endl;
            std::cout << "val is " << std::hex << (unsigned)val << std::endl;
            std::cout << "addr is " << std::hex << addr << std::endl;
            std::cout << "actual val is " << std::hex <<
                (unsigned)mem_val << std::endl;
            return 1;
        }

        return 0;
    }

    static int movw_binary_gen_inddecgen(Sh4 *cpu, Memory *mem,
                                         RandGen32 *randgen32) {
        int failed = 0;

        for (int reg_src = 0; reg_src < 16; reg_src++) {
            for (int reg_dst = 0; reg_dst < 16; reg_dst++) {
                addr32_t addr = randgen32->pick_range(2, mem->get_size() - 2);
                failed = failed ||
                    do_movb_binary_gen_inddecgen(cpu, mem, addr,
                                                 randgen32->pick_val(0),
                                                 reg_src, reg_dst);
            }
        }

        return failed;
    }

    // MOV.L Rm, @-Rn
    // 0010nnnnmmmm0110
    static int do_movl_binary_gen_inddecgen(Sh4 *cpu, Memory *mem,
                                            unsigned addr, uint32_t val,
                                            unsigned reg_src,
                                            unsigned reg_dst) {
        Sh4Prog test_prog;
        std::stringstream ss;
        std::string cmd;
        uint32_t mem_val;

        // increment addr 'cause the opcode is going to decrement it
        addr += 4;

        if (reg_src == reg_dst)
            val = addr;

        ss << "MOV.L R" << reg_src << ", @-R" << reg_dst << "\n";
        cmd = ss.str();
        test_prog.assemble(cmd);
        const Sh4Prog::InstList& inst = test_prog.get_prog();
        mem->load_program(0, inst.begin(), inst.end());

        reset_cpu(cpu);
        *cpu->gen_reg(reg_src) = val;
        *cpu->gen_reg(reg_dst) = addr;
        cpu->exec_inst();

        cpu->read_mem(&mem_val, addr-4, sizeof(mem_val));

        if (reg_src == reg_dst) {
            // special case - val will be decremented because the source and
            // destination are the same register
            val -= 4;
        }

        if (mem_val != val) {
            std::cout << "While running: " << cmd << std::endl;
            std::cout << "val is " << std::hex << (unsigned)val << std::endl;
            std::cout << "addr is " << std::hex << addr << std::endl;
            std::cout << "actual val is " << std::hex <<
                (unsigned)mem_val << std::endl;
            return 1;
        }

        if (*cpu->gen_reg(reg_dst) != addr - 4) {
            std::cout << "While running: " << cmd << std::endl;
            std::cout << "Expected the destination to be decremented "
                "(it was not)" << std::endl;
            std::cout << "val is " << std::hex << (unsigned)val << std::endl;
            std::cout << "addr is " << std::hex << addr << std::endl;
            std::cout << "actual val is " << std::hex <<
                (unsigned)mem_val << std::endl;
            return 1;
        }

        return 0;
    }

    static int movl_binary_gen_inddecgen(Sh4 *cpu, Memory *mem,
                                         RandGen32 *randgen32) {
        int failed = 0;

        for (int reg_src = 0; reg_src < 16; reg_src++) {
            for (int reg_dst = 0; reg_dst < 16; reg_dst++) {
                addr32_t addr = randgen32->pick_range(4, mem->get_size() - 4);
                failed = failed ||
                    do_movb_binary_gen_inddecgen(cpu, mem, addr,
                                                 randgen32->pick_val(0),
                                                 reg_src, reg_dst);
            }
        }

        return failed;
    }

    // MOV.B @Rm+, Rn
    // 0110nnnnmmmm0100
    static int do_movb_binary_indgeninc_gen(Sh4 *cpu, Memory *mem,
                                            unsigned addr, uint8_t val,
                                            unsigned reg_src,
                                            unsigned reg_dst) {
        Sh4Prog test_prog;
        std::stringstream ss;
        std::string cmd;

        if (reg_src == reg_dst)
            val = addr;

        ss << "MOV.B @R" << reg_src << "+, R" << reg_dst << "\n";
        cmd = ss.str();
        test_prog.assemble(cmd);
        const Sh4Prog::InstList& inst = test_prog.get_prog();
        mem->load_program(0, inst.begin(), inst.end());

        reset_cpu(cpu);
        *cpu->gen_reg(reg_src) = addr;
        cpu->write_mem(&val, addr, sizeof(val));
        cpu->exec_inst();

        if (*cpu->gen_reg(reg_dst) != int32_t(val)) {
            std::cout << "While running: " << cmd << std::endl;
            std::cout << "val is " << std::hex << (unsigned)val << std::endl;
            std::cout << "addr is " << std::hex << addr << std::endl;
            std::cout << "actual val is " << std::hex <<
                *cpu->gen_reg(reg_dst) << std::endl;
            return 1;
        }

        if (*cpu->gen_reg(reg_src) != 1 + addr) {
            std::cout << "While running: " << cmd << std::endl;
            std::cout << "The source register did not incrment properly" << std::endl;
            std::cout << "val is " << std::hex << (unsigned)val << std::endl;
            std::cout << "addr is " << std::hex << addr << std::endl;
            std::cout << "actual val is " << std::hex <<
                *cpu->gen_reg(reg_dst) << std::endl;
        }

        return 0;
    }

    static int movb_binary_indgeninc_gen(Sh4 *cpu, Memory *mem,
                                         RandGen32 *randgen32) {
        int failed = 0;

        for (int reg_src = 0; reg_src < 16; reg_src++) {
            for (int reg_dst = 0; reg_dst < 16; reg_dst++) {
                addr32_t addr = randgen32->pick_range(0, mem->get_size() - 2);
                failed = failed ||
                    do_movb_binary_gen_inddecgen(cpu, mem, addr,
                                                 randgen32->pick_val(0),
                                                 reg_src, reg_dst);
            }
        }

        return failed;
    }

    // MOV.W @Rm+, Rn
    // 0110nnnnmmmm0101
    static int do_movw_binary_indgeninc_gen(Sh4 *cpu, Memory *mem,
                                            unsigned addr, uint16_t val,
                                            unsigned reg_src,
                                            unsigned reg_dst) {
        Sh4Prog test_prog;
        std::stringstream ss;
        std::string cmd;

        if (reg_src == reg_dst)
            val = addr;

        ss << "MOV.W @R" << reg_src << "+, R" << reg_dst << "\n";
        cmd = ss.str();
        test_prog.assemble(cmd);
        const Sh4Prog::InstList& inst = test_prog.get_prog();
        mem->load_program(0, inst.begin(), inst.end());

        reset_cpu(cpu);
        *cpu->gen_reg(reg_src) = addr;
        cpu->write_mem(&val, addr, sizeof(val));
        cpu->exec_inst();

        if (*cpu->gen_reg(reg_dst) != int32_t(val)) {
            std::cout << "While running: " << cmd << std::endl;
            std::cout << "val is " << std::hex << (unsigned)val << std::endl;
            std::cout << "addr is " << std::hex << addr << std::endl;
            std::cout << "actual val is " << std::hex <<
                *cpu->gen_reg(reg_dst) << std::endl;
            return 1;
        }

        if (*cpu->gen_reg(reg_src) != 2 + addr) {
            std::cout << "While running: " << cmd << std::endl;
            std::cout << "The source register did not incrment properly" << std::endl;
            std::cout << "val is " << std::hex << (unsigned)val << std::endl;
            std::cout << "addr is " << std::hex << addr << std::endl;
            std::cout << "actual val is " << std::hex <<
                *cpu->gen_reg(reg_dst) << std::endl;
        }

        return 0;
    }

    static int movw_binary_indgeninc_gen(Sh4 *cpu, Memory *mem,
                                         RandGen32 *randgen32) {
        int failed = 0;

        for (int reg_src = 0; reg_src < 16; reg_src++) {
            for (int reg_dst = 0; reg_dst < 16; reg_dst++) {
                addr32_t addr = randgen32->pick_range(0, mem->get_size() - 3);
                failed = failed ||
                    do_movw_binary_gen_inddecgen(cpu, mem,
                                                 randgen32->pick_val(0) %
                                                 mem->get_size(),
                                                 randgen32->pick_val(0),
                                                 reg_src, reg_dst);
            }
        }

        return failed;
    }

    // MOV.L @Rm+, Rn
    // 0110nnnnmmmm0110
    static int do_movl_binary_indgeninc_gen(Sh4 *cpu, Memory *mem,
                                            unsigned addr, uint32_t val,
                                            unsigned reg_src,
                                            unsigned reg_dst) {
        Sh4Prog test_prog;
        std::stringstream ss;
        std::string cmd;

        if (reg_src == reg_dst)
            val = addr;

        ss << "MOV.L @R" << reg_src << "+, R" << reg_dst << "\n";
        cmd = ss.str();
        test_prog.assemble(cmd);
        const Sh4Prog::InstList& inst = test_prog.get_prog();
        mem->load_program(0, inst.begin(), inst.end());

        reset_cpu(cpu);
        *cpu->gen_reg(reg_src) = addr;
        cpu->write_mem(&val, addr, sizeof(val));
        cpu->exec_inst();

        if (*cpu->gen_reg(reg_dst) != int32_t(val)) {
            std::cout << "While running: " << cmd << std::endl;
            std::cout << "val is " << std::hex << (unsigned)val << std::endl;
            std::cout << "addr is " << std::hex << addr << std::endl;
            std::cout << "actual val is " << std::hex <<
                *cpu->gen_reg(reg_dst) << std::endl;
            return 1;
        }

        if (*cpu->gen_reg(reg_src) != 4 + addr) {
            std::cout << "While running: " << cmd << std::endl;
            std::cout << "The source register did not incrment properly" << std::endl;
            std::cout << "val is " << std::hex << (unsigned)val << std::endl;
            std::cout << "addr is " << std::hex << addr << std::endl;
            std::cout << "actual val is " << std::hex <<
                *cpu->gen_reg(reg_dst) << std::endl;
        }

        return 0;
    }

    static int movl_binary_indgeninc_gen(Sh4 *cpu, Memory *mem,
                                         RandGen32 *randgen32) {
        int failed = 0;

        for (int reg_src = 0; reg_src < 16; reg_src++) {
            for (int reg_dst = 0; reg_dst < 16; reg_dst++) {
                addr32_t addr = randgen32->pick_range(0, mem->get_size() - 5);
                failed = failed ||
                    do_movl_binary_gen_inddecgen(cpu, mem, addr,
                                                 randgen32->pick_val(0),
                                                 reg_src, reg_dst);
            }
        }

        return failed;
    }

    // MOV.B R0, @(disp, Rn)
    // 10000000nnnndddd
    static int movb_binary_r0_binind_disp_gen(Sh4 *cpu, Memory *mem,
                                              RandGen32 *randgen32) {
        int failed = 0;

        for (int reg_no = 0; reg_no < 16; reg_no++) {
            for (int disp = 0; disp < 4; disp++) {
                reg32_t base =
                    randgen32->pick_range(0, mem->get_size() - 1 - 0xf);
                DispArg disp_arg(cpu, mem, disp);
                GenRegArg genreg_arg(cpu, mem, reg_no, base);
                uint32_t val = randgen32->pick_val(0);
                if (reg_no == 0) // if both genregs are the same register
                    val = base;
                BinIndArg<DispArg, GenRegArg, 1> dst_op(
                    cpu, mem, disp_arg, genreg_arg,
                    randgen32->pick_val(0), val & 0xff);
                GenRegArg r0_op(cpu, mem, 0, val);
                failed = failed ||
                    do_binary_reg_reg(cpu, mem, "MOV.B", r0_op, dst_op,
                                      DefValidationFunc());
            }
        }

        return failed;
    }

    // MOV.W R0, @(disp, Rn)
    // 10000001nnnndddd
    static int movw_binary_r0_binind_disp_gen(Sh4 *cpu, Memory *mem,
                                              RandGen32 *randgen32) {
        int failed = 0;

        for (int reg_no = 0; reg_no < 16; reg_no++) {
            for (int disp = 0; disp < 4; disp++) {
                reg32_t base =
                    randgen32->pick_range(0, mem->get_size() - 1 - 0xf);
                DispArg disp_arg(cpu, mem, disp);
                GenRegArg genreg_arg(cpu, mem, reg_no, base);
                uint32_t val = randgen32->pick_val(0);
                if (reg_no == 0) // if both genregs are the same register
                    val = base;
                BinIndArg<DispArg, GenRegArg, 2> dst_op(
                    cpu, mem, disp_arg, genreg_arg,
                    randgen32->pick_val(0), val & 0xffff, 2);
                GenRegArg r0_op(cpu, mem, 0, val);
                failed = failed ||
                    do_binary_reg_reg(cpu, mem, "MOV.W", r0_op, dst_op,
                                      DefValidationFunc());
            }
        }

        return failed;
    }

    static int do_movl_binary_gen_binind_disp_gen(Sh4 *cpu, Memory *mem,
                                                  uint8_t disp, reg32_t base,
                                                  uint32_t val, int reg_base,
                                                  int reg_src) {
        std::stringstream ss;
        std::string cmd;
        Sh4Prog test_prog;

        if (reg_base == reg_src) {
            val = base;
        }

        ss << "MOV.L R" << reg_src << ", @(" << (int)disp << ", R" <<
            reg_base << ")\n";
        cmd = ss.str();
        test_prog.assemble(cmd);
        const Sh4Prog::InstList& inst = test_prog.get_prog();
        mem->load_program(0, inst.begin(), inst.end());

        reset_cpu(cpu);
        *cpu->gen_reg(reg_src) = val;
        *cpu->gen_reg(reg_base) = base;
        cpu->exec_inst();

        uint32_t mem_val;
        cpu->read_mem(&mem_val, disp * 4 + base, sizeof(mem_val));
        if (mem_val != val) {
            std::cout << "While running: " << cmd << std::endl;
            std::cout << "val is " << std::hex << (unsigned)val << std::endl;
            std::cout << "disp is " << std::hex << (unsigned)disp << std::endl;
            std::cout << "base is " << std::hex << base << std::endl;
            std::cout << "actual val is " << std::hex << (unsigned)mem_val <<
                std::endl;
            return 1;
        }

        return 0;
    }

    static int movl_binary_gen_binind_disp_gen(Sh4 *cpu, Memory *mem,
                                               RandGen32 *randgen32) {
        int failed = 0;

        for (int reg_src = 0; reg_src < 16; reg_src++) {
            for (int reg_base = 0; reg_base < 16; reg_base++) {
                for (int disp = 0; disp < 4; disp++) {
                    reg32_t base =
                        randgen32->pick_range(0, mem->get_size() - 1 - 0xf);
                    DispArg disp_arg(cpu, mem, disp);
                    GenRegArg genreg_base_arg(cpu, mem, reg_base, base);
                    uint32_t val = randgen32->pick_val(0);
                    if (reg_base == reg_src)
                        val = base; // both genregs are the same register
                    BinIndArg<DispArg, GenRegArg, 4> dst_op(
                        cpu, mem, disp_arg, genreg_base_arg,
                        randgen32->pick_val(0), val, 4);
                    GenRegArg reg_src_op(cpu, mem, reg_src, val);
                    failed = failed ||
                        do_binary_reg_reg(cpu, mem, "MOV.L", reg_src_op, dst_op,
                                          DefValidationFunc());
                }
            }
        }

        return failed;
    }

    // MOV.B @(disp, Rm), R0
    // 10000100mmmmdddd
    static int movb_binary_binind_disp_gen_r0(Sh4 *cpu, Memory *mem,
                                              RandGen32 *randgen32) {
        int failed = 0;

        for (int reg_base = 0; reg_base < 16; reg_base++) {
            for (int disp = 0; disp < 4; disp++) {
                reg32_t base =
                    randgen32->pick_range(0, mem->get_size() - 1 - 0xf);
                reg32_t base_final = base;
                if (reg_base == 0)
                    base_final = int32_t(int8_t(base & 0xff));
                DispArg disp_arg(cpu, mem, disp);
                GenRegArg reg_base_arg(cpu, mem, reg_base, base, base_final);
                uint8_t val = randgen32->pick_val(0) & 0xff;
                uint32_t dstreg_init_val = randgen32->pick_val(0);
                if (reg_base == 0) // if both genregs are the same register
                    dstreg_init_val = base;
                BinIndArg<DispArg, GenRegArg, 1> addr_op(
                    cpu, mem, disp_arg, reg_base_arg,
                    val, val, 1, 1, 0, reg_base == 0);
                GenRegArg r0_op(cpu, mem, 0, dstreg_init_val,
                                int32_t(int8_t(val & 0xff)));
                failed = failed ||
                    do_binary_reg_reg(cpu, mem, "MOV.B", addr_op, r0_op,
                                      DefValidationFunc());
            }
        }
        return failed;
    }

    // MOV.W @(disp, Rm), R0
    // 10000101mmmmdddd
    static int movw_binary_binind_disp_gen_r0(Sh4 *cpu, Memory *mem,
                                              RandGen32 *randgen32) {
        int failed = 0;

        for (int reg_base = 0; reg_base < 16; reg_base++) {
            for (int disp = 0; disp < 4; disp++) {
                reg32_t base =
                    randgen32->pick_range(0, mem->get_size() - 2 - 0xf * 2);
                reg32_t base_final = base;
                if (reg_base == 0)
                    base_final = int32_t(int16_t(base & 0xff));
                DispArg disp_arg(cpu, mem, disp);
                GenRegArg reg_base_arg(cpu, mem, reg_base, base, base_final);
                uint16_t val = randgen32->pick_val(0) & 0xffff;
                uint32_t dstreg_init_val = randgen32->pick_val(0);
                if (reg_base == 0) // if both genregs are the same register
                    dstreg_init_val = base;
                BinIndArg<DispArg, GenRegArg, 2> addr_op(
                    cpu, mem, disp_arg, reg_base_arg,
                    val, val, 2, 1, 0, reg_base == 0);
                GenRegArg r0_op(cpu, mem, 0, dstreg_init_val,
                                int32_t(int16_t(val & 0xffff)));
                failed = failed ||
                    do_binary_reg_reg(cpu, mem, "MOV.W", addr_op, r0_op,
                                      DefValidationFunc());
            }
        }
        return failed;
    }

    // MOV.L @(disp, Rm), Rn
    // 0101nnnnmmmmdddd
    static int movl_binary_binind_disp_gen_gen(Sh4 *cpu, Memory *mem,
                                               RandGen32 *randgen32) {
        int failed = 0;

        for (int reg_base = 0; reg_base < 16; reg_base++) {
            for (int reg_dst = 0; reg_dst < 16; reg_dst++) {
                for (int disp = 0; disp < 4; disp++) {
                    reg32_t base =
                        randgen32->pick_range(0, mem->get_size() - 4 - 0xf * 4);
                    reg32_t base_final = base;
                    DispArg disp_arg(cpu, mem, disp);
                    GenRegArg reg_base_arg(cpu, mem, reg_base,
                                           base, base_final);
                    uint32_t val = randgen32->pick_val(0);
                    uint32_t dstreg_init_val = randgen32->pick_val(0);
                    if (reg_base == reg_dst)
                        dstreg_init_val = base;
                    BinIndArg<DispArg, GenRegArg, 4> addr_op(
                        cpu, mem, disp_arg, reg_base_arg,
                        val, val, 4, 1, 0, reg_base == reg_dst);
                    GenRegArg reg_dst_op(cpu, mem, reg_dst,
                                         dstreg_init_val, val);
                    failed = failed ||
                        do_binary_reg_reg(cpu, mem, "MOV.L", addr_op,
                                          reg_dst_op, DefValidationFunc());
                }
            }
        }
        return failed;
    }

    // MOV.B Rm, @(R0, Rn)
    // 0000nnnnmmmm0100
    static int movb_gen_binind_r0_gen(Sh4 *cpu, Memory *mem,
                                      RandGen32 *randgen32) {
        int failed = 0;

        for (int reg_src = 0; reg_src < 16; reg_src++) {
            for (int reg_dst = 0; reg_dst < 16; reg_dst++) {
                reg32_t r0_addr =
                    randgen32->pick_range(0, mem->get_size() - 1 - 0xf) >> 1;

                reg32_t reg_base_addr =
                    randgen32->pick_range(0, mem->get_size() - 1 - 0xf) >> 1;
                if (reg_dst == 0)
                    reg_base_addr = r0_addr;

                reg32_t src_val = randgen32->pick_val(0);
                if (reg_src == reg_dst)
                    src_val = reg_base_addr;
                else if (reg_src == 0)
                    src_val = r0_addr;
                reg32_t expected_final_dst_val = src_val & 0xff;

                reg32_t initial_dst_val = randgen32->pick_val(0);

                GenRegArg reg_src_arg(cpu, mem, reg_src, src_val);
                GenRegArg reg_addr_right_arg(cpu, mem, reg_dst, reg_base_addr);
                GenRegArg reg_addr_left_arg(cpu, mem, 0, r0_addr);

                BinIndArg<GenRegArg, GenRegArg, 1> addr_op(
                    cpu, mem, reg_addr_left_arg, reg_addr_right_arg,
                    initial_dst_val, expected_final_dst_val, 1, 1, 0, false);

                failed = failed ||
                    do_binary_reg_reg(cpu, mem, "MOV.B", reg_src_arg, addr_op,
                                      DefValidationFunc());
            }
        }
        return failed;
    }

    // MOV.W Rm, @(R0, Rn)
    // 0000nnnnmmmm0101
    static int movw_gen_binind_r0_gen(Sh4 *cpu, Memory *mem,
                                      RandGen32 *randgen32) {
        int failed = 0;

        for (int reg_src = 0; reg_src < 16; reg_src++) {
            for (int reg_dst = 0; reg_dst < 16; reg_dst++) {
                reg32_t r0_addr =
                    randgen32->pick_range(0, mem->get_size() - 1 - 0xf) >> 1;

                reg32_t reg_base_addr =
                    randgen32->pick_range(0, mem->get_size() - 1 - 0xf) >> 1;
                if (reg_dst == 0)
                    reg_base_addr = r0_addr;

                reg32_t src_val = randgen32->pick_val(0);
                if (reg_src == reg_dst)
                    src_val = reg_base_addr;
                else if (reg_src == 0)
                    src_val = r0_addr;
                reg32_t expected_final_dst_val = src_val & 0xffff;

                reg32_t initial_dst_val = randgen32->pick_val(0);

                GenRegArg reg_src_arg(cpu, mem, reg_src, src_val);
                GenRegArg reg_addr_right_arg(cpu, mem, reg_dst, reg_base_addr);
                GenRegArg reg_addr_left_arg(cpu, mem, 0, r0_addr);

                BinIndArg<GenRegArg, GenRegArg, 2> addr_op(
                    cpu, mem, reg_addr_left_arg, reg_addr_right_arg,
                    initial_dst_val, expected_final_dst_val, 1, 1, 0, false);

                failed = failed ||
                    do_binary_reg_reg(cpu, mem, "MOV.W", reg_src_arg, addr_op,
                                      DefValidationFunc());
            }
        }
        return failed;
    }

    // MOV.L Rm, @(R0, Rn)
    // 0001nnnnmmmmdddd
    static int movl_gen_binind_r0_gen(Sh4 *cpu, Memory *mem,
                                      RandGen32 *randgen32) {
        int failed = 0;

        for (int reg_src = 0; reg_src < 16; reg_src++) {
            for (int reg_dst = 0; reg_dst < 16; reg_dst++) {
                reg32_t r0_addr =
                    randgen32->pick_range(0, mem->get_size() - 1 - 0xf) >> 1;

                reg32_t reg_base_addr =
                    randgen32->pick_range(0, mem->get_size() - 1 - 0xf) >> 1;
                if (reg_dst == 0)
                    reg_base_addr = r0_addr;

                reg32_t src_val = randgen32->pick_val(0);
                if (reg_src == reg_dst)
                    src_val = reg_base_addr;
                else if (reg_src == 0)
                    src_val = r0_addr;
                reg32_t expected_final_dst_val = src_val;

                reg32_t initial_dst_val = randgen32->pick_val(0);

                GenRegArg reg_src_arg(cpu, mem, reg_src, src_val);
                GenRegArg reg_addr_right_arg(cpu, mem, reg_dst, reg_base_addr);
                GenRegArg reg_addr_left_arg(cpu, mem, 0, r0_addr);

                BinIndArg<GenRegArg, GenRegArg, 4> addr_op(
                    cpu, mem, reg_addr_left_arg, reg_addr_right_arg,
                    initial_dst_val, expected_final_dst_val, 1, 1, 0, false);

                failed = failed ||
                    do_binary_reg_reg(cpu, mem, "MOV.L", reg_src_arg, addr_op,
                                      DefValidationFunc());
            }
        }
        return failed;
    }

    // MOV.B @(R0, Rm), Rn
    // 0000nnnnmmmm1100
    static int binary_movb_binind_r0_gen_gen(Sh4 *cpu, Memory *mem,
                                             RandGen32 *randgen32) {
        int failed = 0;

        for (int reg_base = 0; reg_base < 16; reg_base++) {
            for (int reg_dst = 0; reg_dst < 16; reg_dst++) {
                reg32_t r0_addr =
                    randgen32->pick_range(0, mem->get_size() - 1 - 0xf) >> 1;

                reg32_t reg_base_addr =
                    randgen32->pick_range(0, mem->get_size() - 1 - 0xf) >> 1;
                if (reg_base == 0)
                    reg_base_addr = r0_addr;

                reg32_t src_val = randgen32->pick_val(0) & 0xff;

                reg32_t initial_dst_val = randgen32->pick_val(0);
                if (reg_dst == 0)
                    initial_dst_val = r0_addr;
                else if (reg_dst == reg_base)
                    initial_dst_val = reg_base_addr;

                reg32_t expected_final_dst_val = int32_t(int8_t(src_val));

                GenRegArg reg_dst_arg(cpu, mem, reg_dst, initial_dst_val,
                                      expected_final_dst_val);
                GenRegArg reg_addr_right_arg(cpu, mem, reg_base, reg_base_addr);
                GenRegArg reg_addr_left_arg(cpu, mem, 0, r0_addr);

                BinIndArg<GenRegArg, GenRegArg, 1> addr_op(
                    cpu, mem, reg_addr_left_arg, reg_addr_right_arg,
                    src_val, src_val, 1, 1, 0,
                    reg_base == 0 || reg_dst == 0 || reg_base == reg_dst);

                failed = failed ||
                    do_binary_reg_reg(cpu, mem, "MOV.B", addr_op, reg_dst_arg,
                                      DefValidationFunc());
            }
        }
        return failed;
    }

    // MOV.W @(R0, Rm), Rn
    // 0000nnnnmmmm1101
    static int binary_movw_binind_r0_gen_gen(Sh4 *cpu, Memory *mem,
                                             RandGen32 *randgen32) {
        int failed = 0;

        for (int reg_base = 0; reg_base < 16; reg_base++) {
            for (int reg_dst = 0; reg_dst < 16; reg_dst++) {
                reg32_t r0_addr =
                    randgen32->pick_range(0, mem->get_size() - 1 - 0xf) >> 1;

                reg32_t reg_base_addr =
                    randgen32->pick_range(0, mem->get_size() - 1 - 0xf) >> 1;
                if (reg_base == 0)
                    reg_base_addr = r0_addr;

                reg32_t src_val = randgen32->pick_val(0) & 0xffff;

                reg32_t initial_dst_val = randgen32->pick_val(0);
                if (reg_dst == 0)
                    initial_dst_val = r0_addr;
                else if (reg_dst == reg_base)
                    initial_dst_val = reg_base_addr;

                reg32_t expected_final_dst_val = int32_t(int16_t(src_val));

                GenRegArg reg_dst_arg(cpu, mem, reg_dst, initial_dst_val,
                                      expected_final_dst_val);
                GenRegArg reg_addr_right_arg(cpu, mem, reg_base, reg_base_addr);
                GenRegArg reg_addr_left_arg(cpu, mem, 0, r0_addr);

                BinIndArg<GenRegArg, GenRegArg, 2> addr_op(
                    cpu, mem, reg_addr_left_arg, reg_addr_right_arg,
                    src_val, src_val, 1, 1, 0,
                    reg_base == 0 || reg_dst == 0 || reg_base == reg_dst);

                failed = failed ||
                    do_binary_reg_reg(cpu, mem, "MOV.W", addr_op, reg_dst_arg,
                                      DefValidationFunc());
            }
        }
        return failed;
    }

    // MOV.L @(R0, Rm), Rn
    // 0000nnnnmmmm1110
    static int binary_movl_binind_r0_gen_gen(Sh4 *cpu, Memory *mem,
                                             RandGen32 *randgen32) {
        int failed = 0;

        for (int reg_base = 0; reg_base < 16; reg_base++) {
            for (int reg_dst = 0; reg_dst < 16; reg_dst++) {
                reg32_t r0_addr =
                    randgen32->pick_range(0, mem->get_size() - 1 - 0xf) >> 1;

                reg32_t reg_base_addr =
                    randgen32->pick_range(0, mem->get_size() - 1 - 0xf) >> 1;
                if (reg_base == 0)
                    reg_base_addr = r0_addr;

                reg32_t src_val = randgen32->pick_val(0);

                reg32_t initial_dst_val = randgen32->pick_val(0);
                if (reg_dst == 0)
                    initial_dst_val = r0_addr;
                else if (reg_dst == reg_base)
                    initial_dst_val = reg_base_addr;

                reg32_t expected_final_dst_val = src_val;

                GenRegArg reg_dst_arg(cpu, mem, reg_dst, initial_dst_val,
                                      expected_final_dst_val);
                GenRegArg reg_addr_right_arg(cpu, mem, reg_base, reg_base_addr);
                GenRegArg reg_addr_left_arg(cpu, mem, 0, r0_addr);

                BinIndArg<GenRegArg, GenRegArg, 4> addr_op(
                    cpu, mem, reg_addr_left_arg, reg_addr_right_arg,
                    src_val, src_val, 1, 1, 0,
                    reg_base == 0 || reg_dst == 0 || reg_base == reg_dst);

                failed = failed ||
                    do_binary_reg_reg(cpu, mem, "MOV.L", addr_op, reg_dst_arg,
                                      DefValidationFunc());
            }
        }
        return failed;
    }

    // MOV.B R0, @(disp, GBR)
    // 11000000dddddddd
    static int binary_movb_r0_binind_disp_gbr(Sh4 *cpu, Memory *mem,
                                              RandGen32 *randgen32) {
        int failed = 0;

        for (unsigned disp = 0; disp < 256; disp++) {
            unsigned gbr_max = mem->get_size() - 1 - disp;
            reg32_t gbr_val = randgen32->pick_range(0, gbr_max);

            uint8_t src_val = randgen32->pick_val(0);

            GenRegArg r0_arg(cpu, mem, 0, src_val);
            DispArg disp_arg(cpu, mem, disp);
            SpecRegArg gbr_arg(cpu, mem, "GBR", &cpu->reg.gbr,  gbr_val);

            BinIndArg<DispArg, SpecRegArg, 1> addr_arg(cpu, mem, disp_arg,
                                                       gbr_arg,
                                                       randgen32->pick_val(0),
                                                       src_val, 1, 1);
            failed = failed ||
                do_binary_reg_reg(cpu, mem, "MOV.B", r0_arg, addr_arg,
                                  DefValidationFunc());
        }

        return failed;
    }

    // MOV.W R0, @(disp, GBR)
    // 11000001dddddddd
    static int binary_movw_r0_binind_disp_gbr(Sh4 *cpu, Memory *mem,
                                              RandGen32 *randgen32) {
        int failed = 0;

        for (unsigned disp = 0; disp < 256; disp++) {
            unsigned gbr_max = mem->get_size() - 2 - disp * 2;
            reg32_t gbr_val = randgen32->pick_range(0, gbr_max);

            uint16_t src_val = randgen32->pick_val(0);

            GenRegArg r0_arg(cpu, mem, 0, src_val);
            DispArg disp_arg(cpu, mem, disp);
            SpecRegArg gbr_arg(cpu, mem, "GBR", &cpu->reg.gbr,  gbr_val);

            BinIndArg<DispArg, SpecRegArg, 2> addr_arg(cpu, mem, disp_arg,
                                                       gbr_arg,
                                                       randgen32->pick_val(0),
                                                       src_val, 2, 1);
            failed = failed ||
                do_binary_reg_reg(cpu, mem, "MOV.W", r0_arg, addr_arg,
                                  DefValidationFunc());
        }

        return failed;
    }

    // MOV.L R0, @(disp, GBR)
    // 11000010dddddddd
    static int binary_movl_r0_binind_disp_gbr(Sh4 *cpu, Memory *mem,
                                              RandGen32 *randgen32) {
        int failed = 0;

        for (unsigned disp = 0; disp < 256; disp++) {
            unsigned gbr_max = mem->get_size() - 4 - disp * 4;
            reg32_t gbr_val = randgen32->pick_range(0, gbr_max);

            uint32_t src_val = randgen32->pick_val(0);

            GenRegArg r0_arg(cpu, mem, 0, src_val);
            DispArg disp_arg(cpu, mem, disp);
            SpecRegArg gbr_arg(cpu, mem, "GBR", &cpu->reg.gbr,  gbr_val);

            BinIndArg<DispArg, SpecRegArg, 4> addr_arg(cpu, mem, disp_arg,
                                                       gbr_arg,
                                                       randgen32->pick_val(0),
                                                       src_val, 4, 1);
            failed = failed ||
                do_binary_reg_reg(cpu, mem, "MOV.L", r0_arg, addr_arg,
                                  DefValidationFunc());
        }

        return failed;
    }

    // MOV.B @(disp, GBR), R0
    // 11000100dddddddd
    static int binary_movb_binind_disp_gbr_r0(Sh4 *cpu, Memory *mem,
                                              RandGen32 *randgen32) {
        int failure = 0;

        for (unsigned disp = 0; disp < 256; disp++) {
            unsigned gbr_max = mem->get_size() - 1 - disp;
            reg32_t gbr_val = randgen32->pick_range(0, gbr_max);

            uint8_t src_val = randgen32->pick_val(0);

            GenRegArg r0_arg(cpu, mem, 0, randgen32->pick_val(0),
                             int32_t(int8_t(src_val)));
            DispArg disp_arg(cpu, mem, disp);
            SpecRegArg gbr_arg(cpu, mem, "GBR", &cpu->reg.gbr, gbr_val);

            BinIndArg<DispArg, SpecRegArg, 1> addr_arg(cpu, mem, disp_arg,
                                                       gbr_arg, src_val,
                                                       src_val, 1, 1);
            failure = failure ||
                do_binary_reg_reg(cpu, mem, "MOV.B", addr_arg, r0_arg,
                                  DefValidationFunc());
        }

        return failure;
    }

    // MOV.W @(disp, GBR), R0
    // 11000101dddddddd
    static int binary_movw_binind_disp_gbr_r0(Sh4 *cpu, Memory *mem,
                                              RandGen32 *randgen32) {
        int failure = 0;

        for (unsigned disp = 0; disp < 256; disp++) {
            unsigned gbr_max = mem->get_size() - 2 - disp * 2;
            reg32_t gbr_val = randgen32->pick_range(0, gbr_max);

            uint8_t src_val = randgen32->pick_val(0);

            GenRegArg r0_arg(cpu, mem, 0, randgen32->pick_val(0),
                             int32_t(int16_t(src_val)));
            DispArg disp_arg(cpu, mem, disp);
            SpecRegArg gbr_arg(cpu, mem, "GBR", &cpu->reg.gbr, gbr_val);

            BinIndArg<DispArg, SpecRegArg, 2> addr_arg(cpu, mem, disp_arg,
                                                       gbr_arg, src_val,
                                                       src_val, 2, 1);
            failure = failure ||
                do_binary_reg_reg(cpu, mem, "MOV.W", addr_arg, r0_arg,
                                  DefValidationFunc());
        }

        return failure;
    }

    // MOV.L @(disp, GBR), R0
    // 11000110dddddddd
    static int binary_movl_binind_disp_gbr_r0(Sh4 *cpu, Memory *mem,
                                              RandGen32 *randgen32) {
        int failure = 0;

        for (unsigned disp = 0; disp < 256; disp++) {
            unsigned gbr_max = mem->get_size() - 4 - disp * 4;
            reg32_t gbr_val = randgen32->pick_range(0, gbr_max);

            uint8_t src_val = randgen32->pick_val(0);

            GenRegArg r0_arg(cpu, mem, 0, randgen32->pick_val(0), src_val);
            DispArg disp_arg(cpu, mem, disp);
            SpecRegArg gbr_arg(cpu, mem, "GBR", &cpu->reg.gbr, gbr_val);

            BinIndArg<DispArg, SpecRegArg, 4> addr_arg(cpu, mem, disp_arg,
                                                       gbr_arg, src_val,
                                                       src_val, 4, 1);
            failure = failure ||
                do_binary_reg_reg(cpu, mem, "MOV.L", addr_arg, r0_arg,
                                  DefValidationFunc());
        }

        return failure;
    }

    // MOVA @(disp, PC), R0
    // 11000111dddddddd
    static int do_binary_mova_binind_disp_pc_r0(Sh4 *cpu, Memory *mem,
                                                uint8_t disp, reg32_t pc_val) {
        std::stringstream ss;
        std::string cmd;
        Sh4Prog test_prog;

        ss << "MOVA @(" << (unsigned)disp << ", PC), R0\n";
        cmd = ss.str();
        test_prog.assemble(cmd);
        const Sh4Prog::InstList& inst = test_prog.get_prog();
        mem->load_program(pc_val, inst.begin(), inst.end());

        reset_cpu(cpu);
        cpu->reg.pc = pc_val;
        cpu->exec_inst();

        reg32_t expected_val = disp * 4 + (pc_val & ~3) + 4;

        if (*cpu->gen_reg(0) != expected_val) {
            std::cout << "ERROR while running \"" << cmd << "\"" << std::endl;
            std::cout << "expected value was " << std::hex <<
                expected_val << std::endl;
            std::cout << "actual value was " << std::hex << *cpu->gen_reg(0) <<
                std::endl;
            std::cout << "PC value was " << std::hex << pc_val << std::endl;
            return 1;
        }

        return 0;
    }

    static int binary_mova_binind_disp_pc_r0(Sh4 *cpu, Memory *mem,
                                             RandGen32 *randgen32) {
        int failure = 0;

        for (int disp = 0; disp <= 0xff; disp++) {
            reg32_t pc_val = randgen32->pick_range(0, (mem->get_size() - 4 - disp * 4) & ~1);
            failure = failure ||
                do_binary_mova_binind_disp_pc_r0(cpu, mem, disp, pc_val);
        }

        return failure;
    }

    static int binary_ldc_test_wrapper(Sh4 *cpu, Memory *mem,
                                       RandGen32 *randgen32,
                                       char const *dst_name, reg32_t *dstp) {
        int failure = 0;
        for (unsigned reg_no = 0; reg_no < 16; reg_no++) {
            reg32_t src_val = randgen32->pick_val(0);
            reg32_t dst_val = randgen32->pick_val(0);

            if (dstp == &cpu->reg.sr) {
                // TODO: it may be better not to set *dstp at all in this case
                dst_val = Sh4::SR_MD_MASK | (cpu->reg.sr & Sh4::SR_RB_MASK);
                src_val = (src_val & ~Sh4::SR_RB_MASK) |
                    (cpu->reg.sr & Sh4::SR_RB_MASK);
            }

            GenRegArg src_reg = GenRegArg(cpu, mem, reg_no, src_val);
            SpecRegArg dst_reg = SpecRegArg(cpu, mem, dst_name,
                                            dstp, dst_val, src_val);
            failure = failure ||
                do_binary_reg_reg<GenRegArg, SpecRegArg, DefValidationFunc>(
                    cpu, mem, "LDC", src_reg, dst_reg, DefValidationFunc());
        }
        return failure;
    }

    // LDC Rm, SR
    // 0100mmmm00001110
    static int binary_ldc_gen_sr(Sh4 *cpu, Memory *mem, RandGen32 *randgen32) {
        return binary_ldc_test_wrapper(cpu, mem, randgen32,
                                       "SR", &cpu->reg.sr);
    }

    // LDC Rm, GBR
    // 0100mmmm00011110
    static int binary_ldc_gen_gbr(Sh4 *cpu, Memory *mem, RandGen32 *randgen32) {
        return binary_ldc_test_wrapper(cpu, mem, randgen32,
                                       "GBR", &cpu->reg.gbr);
    }

    // LDC Rm, VBR
    // 0100mmmm00101110
    static int binary_ldc_gen_vbr(Sh4 *cpu, Memory *mem, RandGen32 *randgen32) {
        return binary_ldc_test_wrapper(cpu, mem, randgen32,
                                       "VBR", &cpu->reg.vbr);
    }

    // LDC Rm, SSR
    // 0100mmmm00111110
    static int binary_ldc_gen_ssr(Sh4 *cpu, Memory *mem, RandGen32 *randgen32) {
        return binary_ldc_test_wrapper(cpu, mem, randgen32,
                                       "SSR", &cpu->reg.ssr);
    }

    // LDC Rm, SPC
    // 0100mmmm01001110
    static int binary_ldc_gen_spc(Sh4 *cpu, Memory *mem, RandGen32 *randgen32) {
        return binary_ldc_test_wrapper(cpu, mem, randgen32,
                                       "SPC", &cpu->reg.spc);
    }

    // LDC Rm, DBR
    // 0100mmmm11111010
    static int binary_ldc_gen_dbr(Sh4 *cpu, Memory *mem, RandGen32 *randgen32) {
        return binary_ldc_test_wrapper(cpu, mem, randgen32,
                                       "DBR", &cpu->reg.dbr);
    }

    // LDC Rm, Rn_BANK
    // 0100mmmm1nnn1110
    static int do_binary_ldc_gen_bank(Sh4 *cpu, Memory *mem,
                                      unsigned reg_no, unsigned bank_reg_no,
                                      reg32_t reg_val) {
        Sh4Prog test_prog;
        std::stringstream ss;
        std::string cmd;

        ss << "LDC R" << reg_no << ", R" << bank_reg_no << "_BANK\n";
        cmd = ss.str();
        test_prog.assemble(cmd);
        const Sh4Prog::InstList& inst = test_prog.get_prog();
        mem->load_program(0, inst.begin(), inst.end());

        reset_cpu(cpu);
        *cpu->gen_reg(reg_no) = reg_val;
        cpu->exec_inst();

        reg32_t bank_reg_val = *cpu->bank_reg(bank_reg_no);

        if (bank_reg_val != reg_val) {
            std::cout << "While running: " << cmd << std::endl;
            std::cout << "reg_val is " << std::hex << reg_val << std::endl;
            std::cout << "actual val is " << std::hex << bank_reg_val <<
                std::endl;
            return 1;
        }

        return 0;
    }

    static int binary_ldc_gen_bank(Sh4 *cpu, Memory *mem,
                                   RandGen32 *randgen32) {
        int failure = 0;

        for (unsigned reg_no = 0; reg_no < 16; reg_no++) {
            for (unsigned bank_reg_no = 0; bank_reg_no < 8; bank_reg_no++) {
                failure = failure ||
                    do_binary_ldc_gen_bank(cpu, mem, reg_no, bank_reg_no,
                                           randgen32->pick_val(0));
            }
        }

        return failure;
    }

    // LDC.L @Rm+, SR
    // 0100mmmm00000111
    static int do_binary_ldcl_indgeninc_sr(Sh4 *cpu, Memory *mem,
                                           unsigned reg_src, addr32_t addr,
                                           uint32_t val) {
        Sh4Prog test_prog;
        std::stringstream ss;
        std::string cmd;

        ss << "LDC.L @R" << reg_src << "+, SR\n";
        cmd = ss.str();
        test_prog.assemble(cmd);
        const Sh4Prog::InstList& inst = test_prog.get_prog();
        mem->load_program(0, inst.begin(), inst.end());

        reset_cpu(cpu);

        *cpu->gen_reg(reg_src) = addr;
        cpu->write_mem(&val, addr, sizeof(val));

        /*
         * Need to restore the original SR because editing SR can cause us to
         * do things that interfere with the test (such as bank-switching).
         */
        reg32_t old_sr = cpu->reg.sr;
        cpu->exec_inst();
        reg32_t new_sr = cpu->reg.sr;
        cpu->reg.sr = old_sr;

        if ((new_sr != val) || (*cpu->gen_reg(reg_src) != 4 + addr)) {
            std::cout << "ERROR while running " << cmd << std::endl;
            std::cout << "address is " << std::hex << addr << std::endl;
            std::cout << "expected value is " << val << std::endl;
            std::cout << "actual value is " << new_sr << std::endl;
            std::cout << "expected output address is " << (4 + addr) <<
                std::endl;
            std::cout << "actual output address is " <<
                *cpu->gen_reg(reg_src) << std::endl;
            return 1;
        }

        return 0;
    }

    static int binary_ldcl_indgeninc_sr(Sh4 *cpu, Memory *mem,
                                        RandGen32 *randgen32) {
        int failure = 0;

        for (unsigned reg_src = 0; reg_src < 16; reg_src++) {
            addr32_t addr = randgen32->pick_range(0, mem->get_size() - 5);
            addr32_t val = randgen32->pick_val(0);

            failure = failure ||
                do_binary_ldcl_indgeninc_sr(cpu, mem, reg_src, addr, val);
        }

        return failure;
    }

    // LDC.L @Rm+, GBR
    // 0100mmmm00010111
    static int do_binary_ldcl_indgeninc_gbr(Sh4 *cpu, Memory *mem,
                                            unsigned reg_src, addr32_t addr,
                                            uint32_t val) {
        Sh4Prog test_prog;
        std::stringstream ss;
        std::string cmd;

        ss << "LDC.L @R" << reg_src << "+, GBR\n";
        cmd = ss.str();
        test_prog.assemble(cmd);
        const Sh4Prog::InstList& inst = test_prog.get_prog();
        mem->load_program(0, inst.begin(), inst.end());

        reset_cpu(cpu);

        *cpu->gen_reg(reg_src) = addr;
        cpu->write_mem(&val, addr, sizeof(val));

        cpu->exec_inst();

        if ((cpu->reg.gbr != val) || (*cpu->gen_reg(reg_src) != 4 + addr)) {
            std::cout << "ERROR while running " << cmd << std::endl;
            std::cout << "address is " << std::hex << addr << std::endl;
            std::cout << "expected value is " << val << std::endl;
            std::cout << "actual value is " << cpu->reg.gbr << std::endl;
            std::cout << "expected output address is " << (4 + addr) <<
                std::endl;
            std::cout << "actual output address is " <<
                *cpu->gen_reg(reg_src) << std::endl;
            return 1;
        }

        return 0;
    }

    static int binary_ldcl_indgeninc_gbr(Sh4 *cpu, Memory *mem,
                                         RandGen32 *randgen32) {
        int failure = 0;

        for (unsigned reg_src = 0; reg_src < 16; reg_src++) {
            addr32_t addr = randgen32->pick_range(0, mem->get_size() - 5);
            addr32_t val = randgen32->pick_val(0);

            failure = failure ||
                do_binary_ldcl_indgeninc_gbr(cpu, mem, reg_src, addr, val);
        }

        return failure;
    }

    // LDC.L @Rm+, VBR
    // 0100mmmm00100111
    static int do_binary_ldcl_indgeninc_vbr(Sh4 *cpu, Memory *mem,
                                            unsigned reg_src, addr32_t addr,
                                            uint32_t val) {
        Sh4Prog test_prog;
        std::stringstream ss;
        std::string cmd;

        ss << "LDC.L @R" << reg_src << "+, VBR\n";
        cmd = ss.str();
        test_prog.assemble(cmd);
        const Sh4Prog::InstList& inst = test_prog.get_prog();
        mem->load_program(0, inst.begin(), inst.end());

        reset_cpu(cpu);

        *cpu->gen_reg(reg_src) = addr;
        cpu->write_mem(&val, addr, sizeof(val));

        cpu->exec_inst();

        if ((cpu->reg.vbr != val) || (*cpu->gen_reg(reg_src) != 4 + addr)) {
            std::cout << "ERROR while running " << cmd << std::endl;
            std::cout << "address is " << std::hex << addr << std::endl;
            std::cout << "expected value is " << val << std::endl;
            std::cout << "actual value is " << cpu->reg.vbr << std::endl;
            std::cout << "expected output address is " << (4 + addr) <<
                std::endl;
            std::cout << "actual output address is " <<
                *cpu->gen_reg(reg_src) << std::endl;
            return 1;
        }

        return 0;
    }

    static int binary_ldcl_indgeninc_vbr(Sh4 *cpu, Memory *mem,
                                         RandGen32 *randgen32) {
        int failure = 0;

        for (unsigned reg_src = 0; reg_src < 16; reg_src++) {
            addr32_t addr = randgen32->pick_range(0, mem->get_size() - 5);
            addr32_t val = randgen32->pick_val(0);

            failure = failure ||
                do_binary_ldcl_indgeninc_vbr(cpu, mem, reg_src, addr, val);
        }

        return failure;
    }

    // LDC.L @Rm+, SSR
    // 0100mmmm00110111
    static int do_binary_ldcl_indgeninc_ssr(Sh4 *cpu, Memory *mem,
                                            unsigned reg_src, addr32_t addr,
                                            uint32_t val) {
        Sh4Prog test_prog;
        std::stringstream ss;
        std::string cmd;

        ss << "LDC.L @R" << reg_src << "+, SSR\n";
        cmd = ss.str();
        test_prog.assemble(cmd);
        const Sh4Prog::InstList& inst = test_prog.get_prog();
        mem->load_program(0, inst.begin(), inst.end());

        reset_cpu(cpu);

        *cpu->gen_reg(reg_src) = addr;
        cpu->write_mem(&val, addr, sizeof(val));

        cpu->exec_inst();

        if ((cpu->reg.ssr != val) || (*cpu->gen_reg(reg_src) != 4 + addr)) {
            std::cout << "ERROR while running " << cmd << std::endl;
            std::cout << "address is " << std::hex << addr << std::endl;
            std::cout << "expected value is " << val << std::endl;
            std::cout << "actual value is " << cpu->reg.ssr << std::endl;
            std::cout << "expected output address is " << (4 + addr) <<
                std::endl;
            std::cout << "actual output address is " <<
                *cpu->gen_reg(reg_src) << std::endl;
            return 1;
        }

        return 0;
    }

    static int binary_ldcl_indgeninc_ssr(Sh4 *cpu, Memory *mem,
                                         RandGen32 *randgen32) {
        int failure = 0;

        for (unsigned reg_src = 0; reg_src < 16; reg_src++) {
            addr32_t addr = randgen32->pick_range(0, mem->get_size() - 5);
            addr32_t val = randgen32->pick_val(0);

            failure = failure ||
                do_binary_ldcl_indgeninc_ssr(cpu, mem, reg_src, addr, val);
        }

        return failure;
    }

    // LDC.L @Rm+, SPC
    // 0100mmmm01000111
    static int do_binary_ldcl_indgeninc_spc(Sh4 *cpu, Memory *mem,
                                            unsigned reg_src, addr32_t addr,
                                            uint32_t val) {
        Sh4Prog test_prog;
        std::stringstream ss;
        std::string cmd;

        ss << "LDC.L @R" << reg_src << "+, SPC\n";
        cmd = ss.str();
        test_prog.assemble(cmd);
        const Sh4Prog::InstList& inst = test_prog.get_prog();
        mem->load_program(0, inst.begin(), inst.end());

        reset_cpu(cpu);

        *cpu->gen_reg(reg_src) = addr;
        cpu->write_mem(&val, addr, sizeof(val));

        cpu->exec_inst();

        if ((cpu->reg.spc != val) || (*cpu->gen_reg(reg_src) != 4 + addr)) {
            std::cout << "ERROR while running " << cmd << std::endl;
            std::cout << "address is " << std::hex << addr << std::endl;
            std::cout << "expected value is " << val << std::endl;
            std::cout << "actual value is " << cpu->reg.spc << std::endl;
            std::cout << "expected output address is " << (4 + addr) <<
                std::endl;
            std::cout << "actual output address is " <<
                *cpu->gen_reg(reg_src) << std::endl;
            return 1;
        }

        return 0;
    }

    static int binary_ldcl_indgeninc_spc(Sh4 *cpu, Memory *mem,
                                         RandGen32 *randgen32) {
        int failure = 0;

        for (unsigned reg_src = 0; reg_src < 16; reg_src++) {
            addr32_t addr = randgen32->pick_range(0, mem->get_size() - 5);
            addr32_t val = randgen32->pick_val(0);

            failure = failure ||
                do_binary_ldcl_indgeninc_spc(cpu, mem, reg_src, addr, val);
        }

        return failure;
    }

    // LDC.L @Rm+, DBR
    // 0100mmmm11110110
    static int do_binary_ldcl_indgeninc_dbr(Sh4 *cpu, Memory *mem,
                                            unsigned reg_src, addr32_t addr,
                                            uint32_t val) {
        Sh4Prog test_prog;
        std::stringstream ss;
        std::string cmd;

        ss << "LDC.L @R" << reg_src << "+, DBR\n";
        cmd = ss.str();
        test_prog.assemble(cmd);
        const Sh4Prog::InstList& inst = test_prog.get_prog();
        mem->load_program(0, inst.begin(), inst.end());

        reset_cpu(cpu);

        *cpu->gen_reg(reg_src) = addr;
        cpu->write_mem(&val, addr, sizeof(val));

        cpu->exec_inst();

        if ((cpu->reg.dbr != val) || (*cpu->gen_reg(reg_src) != 4 + addr)) {
            std::cout << "ERROR while running " << cmd << std::endl;
            std::cout << "address is " << std::hex << addr << std::endl;
            std::cout << "expected value is " << val << std::endl;
            std::cout << "actual value is " << cpu->reg.dbr << std::endl;
            std::cout << "expected output address is " << (4 + addr) <<
                std::endl;
            std::cout << "actual output address is " <<
                *cpu->gen_reg(reg_src) << std::endl;
            return 1;
        }

        return 0;
    }

    static int binary_ldcl_indgeninc_dbr(Sh4 *cpu, Memory *mem,
                                         RandGen32 *randgen32) {
        int failure = 0;

        for (unsigned reg_src = 0; reg_src < 16; reg_src++) {
            addr32_t addr = randgen32->pick_range(0, mem->get_size() - 5);
            addr32_t val = randgen32->pick_val(0);

            failure = failure ||
                do_binary_ldcl_indgeninc_dbr(cpu, mem, reg_src, addr, val);
        }

        return failure;
    }

    static int binary_stc_test_wrapper(Sh4 *cpu, Memory *mem,
                                       RandGen32 *randgen32,
                                       char const *src_name, reg32_t *srcp) {
        int failure = 0;
        for (unsigned reg_no = 0; reg_no < 16; reg_no++) {
            reg32_t src_val = randgen32->pick_val(0);
            reg32_t dst_val = randgen32->pick_val(0);

            // don't accidentally disable privileged mode or switch banks
            if (srcp == &cpu->reg.sr) {
                src_val |= Sh4::SR_MD_MASK;
                src_val &= ~Sh4::SR_RB_MASK;
            }

            SpecRegArg src_reg = SpecRegArg(cpu, mem, src_name,
                                            srcp, src_val);
            GenRegArg dst_reg = GenRegArg(cpu, mem, reg_no, dst_val, src_val);
            failure = failure ||
                do_binary_reg_reg<SpecRegArg, GenRegArg, DefValidationFunc>(
                    cpu, mem, "STC", src_reg, dst_reg, DefValidationFunc());
        }
        return failure;
    }

    // STC SR, Rn
    // 0000nnnn00000010
    static int binary_stc_sr_gen(Sh4 *cpu, Memory *mem, RandGen32 *randgen32) {
        return binary_stc_test_wrapper(cpu, mem, randgen32,
                                       "SR", &cpu->reg.sr);
    }

    // STC GBR, Rn
    // 0000nnnn00010010
    static int binary_stc_gbr_gen(Sh4 *cpu, Memory *mem, RandGen32 *randgen32) {
        return binary_stc_test_wrapper(cpu, mem, randgen32,
                                       "GBR", &cpu->reg.gbr);
    }

    // STC VBR, Rn
    // 0000nnnn00100010
    static int binary_stc_vbr_gen(Sh4 *cpu, Memory *mem, RandGen32 *randgen32) {
        return binary_stc_test_wrapper(cpu, mem, randgen32,
                                       "VBR", &cpu->reg.vbr);
    }

    // STC SSR, Rn
    // 0000nnnn00110010
    static int binary_stc_ssr_gen(Sh4 *cpu, Memory *mem, RandGen32 *randgen32) {
        return binary_stc_test_wrapper(cpu, mem, randgen32,
                                       "SSR", &cpu->reg.ssr);
    }

    // STC SPC, Rn
    // 0000nnnn01000010
    static int binary_stc_spc_gen(Sh4 *cpu, Memory *mem, RandGen32 *randgen32) {
        return binary_stc_test_wrapper(cpu, mem, randgen32,
                                       "SPC", &cpu->reg.spc);
    }

    // STC SGR, Rn
    // 0000nnnn00111010
    static int binary_stc_sgr_gen(Sh4 *cpu, Memory *mem, RandGen32 *randgen32) {
        return binary_stc_test_wrapper(cpu, mem, randgen32,
                                       "SGR", &cpu->reg.sgr);
    }

    // STC DBR, Rn
    // 0000nnnn11111010
    static int binary_stc_dbr_gen(Sh4 *cpu, Memory *mem, RandGen32 *randgen32) {
        return binary_stc_test_wrapper(cpu, mem, randgen32,
                                       "DBR", &cpu->reg.dbr);
    }

    // STC.L SR, @-Rn
    // 0100nnnn00000011
    static int do_binary_stcl_sr_inddecgen(Sh4 *cpu, Memory *mem,
                                           unsigned reg_no, reg32_t sr_val,
                                           addr32_t addr) {
        Sh4Prog test_prog;
        std::stringstream ss;
        std::string cmd;

        // obviously this needs to run in privileged mode
        sr_val |= Sh4::SR_MD_MASK;

        ss << "STC.L SR, @-R" << reg_no << "\n";
        cmd = ss.str();
        test_prog.assemble(cmd);
        const Sh4Prog::InstList& inst = test_prog.get_prog();
        mem->load_program(0, inst.begin(), inst.end());

        reset_cpu(cpu);

        cpu->reg.sr = sr_val;
        *cpu->gen_reg(reg_no) = addr;
        cpu->exec_inst();

        uint32_t mem_val;
        cpu->read_mem(&mem_val, addr - 4, sizeof(mem_val));

        if (sr_val != mem_val || *cpu->gen_reg(reg_no) != (addr - 4)) {
            std::cout << "ERROR while running " << cmd << std::endl;
            std::cout << "Expected value was " << std::hex << sr_val <<
                std::endl;
            std::cout << "Actual value is " << mem_val << std::endl;
            std::cout << "expected output addr is " << (addr - 4) << std::endl;
            std::cout << "actual output addr is " << *cpu->gen_reg(reg_no) <<
                std::endl;
            return 1;
        }

        return 0;
    }

    static int binary_stcl_sr_inddecgen(Sh4 *cpu, Memory *mem,
                                        RandGen32 *randgen32) {
        int failure = 0;

        for (unsigned reg_no = 0; reg_no < 16; reg_no++) {
            addr32_t addr = randgen32->pick_range(4, mem->get_size() - 4);
            failure = failure ||
                do_binary_stcl_sr_inddecgen(cpu, mem, reg_no,
                                            randgen32->pick_val(0), addr);
        }

        return failure;
    }

    // STC.L GBR, @-Rn
    // 0100nnnn00010011
    static int do_binary_stcl_gbr_inddecgen(Sh4 *cpu, Memory *mem,
                                            unsigned reg_no, reg32_t gbr_val,
                                            addr32_t addr) {
        Sh4Prog test_prog;
        std::stringstream ss;
        std::string cmd;

        ss << "STC.L GBR, @-R" << reg_no << "\n";
        cmd = ss.str();
        test_prog.assemble(cmd);
        const Sh4Prog::InstList& inst = test_prog.get_prog();
        mem->load_program(0, inst.begin(), inst.end());

        reset_cpu(cpu);

        cpu->reg.gbr = gbr_val;
        *cpu->gen_reg(reg_no) = addr;
        cpu->exec_inst();

        uint32_t mem_val;
        cpu->read_mem(&mem_val, addr - 4, sizeof(mem_val));

        if (gbr_val != mem_val || *cpu->gen_reg(reg_no) != (addr - 4)) {
            std::cout << "ERROR while running " << cmd << std::endl;
            std::cout << "Expected value was " << std::hex << gbr_val <<
                std::endl;
            std::cout << "Actual value is " << mem_val << std::endl;
            std::cout << "expected output addr is " << (addr - 4) << std::endl;
            std::cout << "actual output addr is " << *cpu->gen_reg(reg_no) <<
                std::endl;
            return 1;
        }

        return 0;
    }

    static int binary_stcl_gbr_inddecgen(Sh4 *cpu, Memory *mem,
                                         RandGen32 *randgen32) {
        int failure = 0;

        for (unsigned reg_no = 0; reg_no < 16; reg_no++) {
            addr32_t addr = randgen32->pick_range(4, mem->get_size() - 4);
            failure = failure ||
                do_binary_stcl_gbr_inddecgen(cpu, mem, reg_no,
                                             randgen32->pick_val(0), addr);
        }

        return failure;
    }

    // STC.L VBR, @-Rn
    // 01n00nnnn00100011
    static int do_binary_stcl_vbr_inddecgen(Sh4 *cpu, Memory *mem,
                                            unsigned reg_no, reg32_t vbr_val,
                                            addr32_t addr) {
        Sh4Prog test_prog;
        std::stringstream ss;
        std::string cmd;

        ss << "STC.L VBR, @-R" << reg_no << "\n";
        cmd = ss.str();
        test_prog.assemble(cmd);
        const Sh4Prog::InstList& inst = test_prog.get_prog();
        mem->load_program(0, inst.begin(), inst.end());

        reset_cpu(cpu);

        cpu->reg.vbr = vbr_val;
        *cpu->gen_reg(reg_no) = addr;
        cpu->exec_inst();

        uint32_t mem_val;
        cpu->read_mem(&mem_val, addr - 4, sizeof(mem_val));

        if (vbr_val != mem_val || *cpu->gen_reg(reg_no) != (addr - 4)) {
            std::cout << "ERROR while running " << cmd << std::endl;
            std::cout << "Expected value was " << std::hex << vbr_val <<
                std::endl;
            std::cout << "Actual value is " << mem_val << std::endl;
            std::cout << "expected output addr is " << (addr - 4) << std::endl;
            std::cout << "actual output addr is " << *cpu->gen_reg(reg_no) <<
                std::endl;
            return 1;
        }

        return 0;
    }

    static int binary_stcl_vbr_inddecgen(Sh4 *cpu, Memory *mem,
                                         RandGen32 *randgen32) {
        int failure = 0;

        for (unsigned reg_no = 0; reg_no < 16; reg_no++) {
            addr32_t addr = randgen32->pick_range(4, mem->get_size() - 4);
            failure = failure ||
                do_binary_stcl_vbr_inddecgen(cpu, mem, reg_no,
                                             randgen32->pick_val(0), addr);
        }

        return failure;
    }

    // STC.L SSR, @-Rn
    // 0100nnnn00110011
    static int do_binary_stcl_ssr_inddecgen(Sh4 *cpu, Memory *mem,
                                            unsigned reg_no, reg32_t ssr_val,
                                            addr32_t addr) {
        Sh4Prog test_prog;
        std::stringstream ss;
        std::string cmd;

        ss << "STC.L SSR, @-R" << reg_no << "\n";
        cmd = ss.str();
        test_prog.assemble(cmd);
        const Sh4Prog::InstList& inst = test_prog.get_prog();
        mem->load_program(0, inst.begin(), inst.end());

        reset_cpu(cpu);

        cpu->reg.ssr = ssr_val;
        *cpu->gen_reg(reg_no) = addr;
        cpu->exec_inst();

        uint32_t mem_val;
        cpu->read_mem(&mem_val, addr - 4, sizeof(mem_val));

        if (ssr_val != mem_val || *cpu->gen_reg(reg_no) != (addr - 4)) {
            std::cout << "ERROR while running " << cmd << std::endl;
            std::cout << "Expected value was " << std::hex << ssr_val <<
                std::endl;
            std::cout << "Actual value is " << mem_val << std::endl;
            std::cout << "expected output addr is " << (addr - 4) << std::endl;
            std::cout << "actual output addr is " << *cpu->gen_reg(reg_no) <<
                std::endl;
            return 1;
        }

        return 0;
    }

    static int binary_stcl_ssr_inddecgen(Sh4 *cpu, Memory *mem,
                                         RandGen32 *randgen32) {
        int failure = 0;

        for (unsigned reg_no = 0; reg_no < 16; reg_no++) {
            addr32_t addr = randgen32->pick_range(4, mem->get_size() - 4);
            failure = failure ||
                do_binary_stcl_ssr_inddecgen(cpu, mem, reg_no,
                                             randgen32->pick_val(0), addr);
        }

        return failure;
    }

    // STC.L SPC, @-Rn
    // 0100nnnn01000011
    static int do_binary_stcl_spc_inddecgen(Sh4 *cpu, Memory *mem,
                                            unsigned reg_no, reg32_t spc_val,
                                            addr32_t addr) {
        Sh4Prog test_prog;
        std::stringstream ss;
        std::string cmd;

        ss << "STC.L SPC, @-R" << reg_no << "\n";
        cmd = ss.str();
        test_prog.assemble(cmd);
        const Sh4Prog::InstList& inst = test_prog.get_prog();
        mem->load_program(0, inst.begin(), inst.end());

        reset_cpu(cpu);

        cpu->reg.spc = spc_val;
        *cpu->gen_reg(reg_no) = addr;
        cpu->exec_inst();

        uint32_t mem_val;
        cpu->read_mem(&mem_val, addr - 4, sizeof(mem_val));

        if (spc_val != mem_val || *cpu->gen_reg(reg_no) != (addr - 4)) {
            std::cout << "ERROR while running " << cmd << std::endl;
            std::cout << "Expected value was " << std::hex << spc_val <<
                std::endl;
            std::cout << "Actual value is " << mem_val << std::endl;
            std::cout << "expected output addr is " << (addr - 4) << std::endl;
            std::cout << "actual output addr is " << *cpu->gen_reg(reg_no) <<
                std::endl;
            return 1;
        }

        return 0;
    }

    static int binary_stcl_spc_inddecgen(Sh4 *cpu, Memory *mem,
                                         RandGen32 *randgen32) {
        int failure = 0;

        for (unsigned reg_no = 0; reg_no < 16; reg_no++) {
            addr32_t addr = randgen32->pick_range(4, mem->get_size() - 4);
            failure = failure ||
                do_binary_stcl_spc_inddecgen(cpu, mem, reg_no,
                                             randgen32->pick_val(0), addr);
        }

        return failure;
    }

    // STC.L SGR, @-Rn
    // 0100nnnn00110010
    static int do_binary_stcl_sgr_inddecgen(Sh4 *cpu, Memory *mem,
                                            unsigned reg_no, reg32_t sgr_val,
                                            addr32_t addr) {
        Sh4Prog test_prog;
        std::stringstream ss;
        std::string cmd;

        ss << "STC.L SGR, @-R" << reg_no << "\n";
        cmd = ss.str();
        test_prog.assemble(cmd);
        const Sh4Prog::InstList& inst = test_prog.get_prog();
        mem->load_program(0, inst.begin(), inst.end());

        reset_cpu(cpu);

        cpu->reg.sgr = sgr_val;
        *cpu->gen_reg(reg_no) = addr;
        cpu->exec_inst();

        uint32_t mem_val;
        cpu->read_mem(&mem_val, addr - 4, sizeof(mem_val));

        if (sgr_val != mem_val || *cpu->gen_reg(reg_no) != (addr - 4)) {
            std::cout << "ERROR while running " << cmd << std::endl;
            std::cout << "Expected value was " << std::hex << sgr_val <<
                std::endl;
            std::cout << "Actual value is " << mem_val << std::endl;
            std::cout << "expected output addr is " << (addr - 4) << std::endl;
            std::cout << "actual output addr is " << *cpu->gen_reg(reg_no) <<
                std::endl;
            return 1;
        }

        return 0;
    }

    static int binary_stcl_sgr_inddecgen(Sh4 *cpu, Memory *mem,
                                         RandGen32 *randgen32) {
        int failure = 0;

        for (unsigned reg_no = 0; reg_no < 16; reg_no++) {
            addr32_t addr = randgen32->pick_range(4, mem->get_size() - 4);
            failure = failure ||
                do_binary_stcl_sgr_inddecgen(cpu, mem, reg_no,
                                             randgen32->pick_val(0), addr);
        }

        return failure;
    }

    // STC.L DBR, @-Rn
    // 0100nnnn11110010
    static int do_binary_stcl_dbr_inddecgen(Sh4 *cpu, Memory *mem,
                                            unsigned reg_no, reg32_t dbr_val,
                                            addr32_t addr) {
        Sh4Prog test_prog;
        std::stringstream ss;
        std::string cmd;

        ss << "STC.L DBR, @-R" << reg_no << "\n";
        cmd = ss.str();
        test_prog.assemble(cmd);
        const Sh4Prog::InstList& inst = test_prog.get_prog();
        mem->load_program(0, inst.begin(), inst.end());

        reset_cpu(cpu);

        cpu->reg.dbr = dbr_val;
        *cpu->gen_reg(reg_no) = addr;
        cpu->exec_inst();

        uint32_t mem_val;
        cpu->read_mem(&mem_val, addr - 4, sizeof(mem_val));

        if (dbr_val != mem_val || *cpu->gen_reg(reg_no) != (addr - 4)) {
            std::cout << "ERROR while running " << cmd << std::endl;
            std::cout << "Expected value was " << std::hex << dbr_val <<
                std::endl;
            std::cout << "Actual value is " << mem_val << std::endl;
            std::cout << "expected output addr is " << (addr - 4) << std::endl;
            std::cout << "actual output addr is " << *cpu->gen_reg(reg_no) <<
                std::endl;
            return 1;
        }

        return 0;
    }

    static int binary_stcl_dbr_inddecgen(Sh4 *cpu, Memory *mem,
                                         RandGen32 *randgen32) {
        int failure = 0;

        for (unsigned reg_no = 0; reg_no < 16; reg_no++) {
            addr32_t addr = randgen32->pick_range(4, mem->get_size() - 4);
            failure = failure ||
                do_binary_stcl_dbr_inddecgen(cpu, mem, reg_no,
                                             randgen32->pick_val(0), addr);
        }

        return failure;
    }

    // LDC.L @Rm+, Rn_BANK
    // 0100mmmm1nnn0111
    static int do_binary_ldcl_indgeninc_bank(Sh4 *cpu, Memory *mem,
                                             unsigned reg_no,
                                             unsigned bank_reg_no,
                                             addr32_t addr, uint32_t val) {
        Sh4Prog test_prog;
        std::stringstream ss;
        std::string cmd;

        ss << "LDC.L @R" << reg_no << "+, R" << bank_reg_no << "_BANK\n";
        cmd = ss.str();
        test_prog.assemble(cmd);
        const Sh4Prog::InstList& inst = test_prog.get_prog();
        mem->load_program(0, inst.begin(), inst.end());

        reset_cpu(cpu);
        *cpu->gen_reg(reg_no) = addr;
        cpu->write_mem(&val, addr, sizeof(val));
        cpu->exec_inst();

        reg32_t bank_reg_val = *cpu->bank_reg(bank_reg_no);

        if (bank_reg_val != val || *cpu->gen_reg(reg_no) != (addr + 4)) {
            std::cout << "While running: " << cmd << std::endl;
            std::cout << "input address is " << addr << std::endl;
            std::cout << "val is " << std::hex << val << std::endl;
            std::cout << "actual val is " << std::hex << bank_reg_val <<
                std::endl;
            std::cout << "expected output address is " << addr + 4 << std::endl;
            std::cout << "actual output address is " << *cpu->gen_reg(reg_no) <<
                std::endl;
            return 1;
        }

        return 0;
    }

    static int binary_ldcl_indgeninc_bank(Sh4 *cpu, Memory *mem,
                                          RandGen32 *randgen32) {
        int failure = 0;

        for (unsigned reg_no = 0; reg_no < 16; reg_no++) {
            for (unsigned bank_reg_no = 0; bank_reg_no < 8; bank_reg_no++) {
                addr32_t addr = randgen32->pick_range(0, mem->get_size() - 5);
                failure = failure ||
                    do_binary_ldcl_indgeninc_bank(cpu, mem, reg_no, bank_reg_no,
                                                  addr, randgen32->pick_val(0));
            }
        }

        return failure;
    }

    // STC Rm_BANK, Rn
    // 0000nnnn1mmm0010
    static int do_binary_stc_bank_gen(Sh4 *cpu, Memory *mem,
                                      unsigned bank_reg_no, unsigned reg_no,
                                      reg32_t val) {
        Sh4Prog test_prog;
        std::stringstream ss;
        std::string cmd;

        ss << "STC R" << bank_reg_no << "_BANK, R" << reg_no << "\n";
        cmd = ss.str();
        test_prog.assemble(cmd);
        const Sh4Prog::InstList& inst = test_prog.get_prog();
        mem->load_program(0, inst.begin(), inst.end());

        reset_cpu(cpu);

        *cpu->bank_reg(bank_reg_no) = val;
        cpu->exec_inst();

        if (*cpu->gen_reg(reg_no) != val) {
            std::cout << "ERROR while running " << cmd << std::endl;
            std::cout << "Expected value was " << std::hex << val <<
                std::endl;
            std::cout << "Actual value is " << *cpu->gen_reg(reg_no) <<
                std::endl;
            return 1;
        }

        return 0;
    }

    static int binary_stc_bank_gen(Sh4 *cpu, Memory *mem,
                                   RandGen32 *randgen32) {
        int failure = 0;
        for (unsigned reg_no = 0; reg_no < 16; reg_no++) {
            for (unsigned bank_reg_no = 0; bank_reg_no < 8; bank_reg_no++) {
                failure = failure ||
                    do_binary_stc_bank_gen(cpu, mem, bank_reg_no, reg_no,
                                           randgen32->pick_val(0));
            }
        }
        return failure;
    }

    // STC.L Rm_BANK, @-Rn
    // 0100nnnn1mmm0011
    static int do_binary_stcl_bank_inddecgen(Sh4 *cpu, Memory *mem,
                                             unsigned bank_reg_no,
                                             unsigned reg_no, reg32_t val,
                                             addr32_t addr) {
        Sh4Prog test_prog;
        std::stringstream ss;
        std::string cmd;

        ss << "STC.L R" << bank_reg_no << "_BANK, @-R" << reg_no << "\n";
        cmd = ss.str();
        test_prog.assemble(cmd);
        const Sh4Prog::InstList& inst = test_prog.get_prog();
        mem->load_program(0, inst.begin(), inst.end());

        reset_cpu(cpu);

        *cpu->bank_reg(bank_reg_no) = val;
        *cpu->gen_reg(reg_no) = addr;
        cpu->exec_inst();

        uint32_t mem_val;
        cpu->read_mem(&mem_val, addr - 4, sizeof(mem_val));

        if (val != mem_val || *cpu->gen_reg(reg_no) != (addr - 4)) {
            std::cout << "ERROR while running " << cmd << std::endl;
            std::cout << "Expected value was " << std::hex << val <<
                std::endl;
            std::cout << "Actual value is " << mem_val << std::endl;
            std::cout << "expected output addr is " << (addr - 4) << std::endl;
            std::cout << "actual output addr is " << *cpu->gen_reg(reg_no) <<
                std::endl;
            return 1;
        }

        return 0;
    }

    static int binary_stcl_bank_inddecgen(Sh4 *cpu, Memory *mem,
                                          RandGen32 *randgen32) {
        int failure = 0;

        for (unsigned reg_no = 0; reg_no < 16; reg_no++) {
            for (unsigned bank_reg_no = 0; bank_reg_no < 8; bank_reg_no++) {
                addr32_t addr = randgen32->pick_range(4, mem->get_size() - 4);
                failure = failure ||
                    do_binary_stcl_bank_inddecgen(cpu, mem, bank_reg_no, reg_no,
                                                  randgen32->pick_val(0), addr);
            }
        }

        return failure;
    }

    static int binary_lds_test_wrapper(Sh4 *cpu, Memory *mem,
                                       RandGen32 *randgen32,
                                       char const *dst_name, reg32_t *dstp) {
        int failure = 0;
        for (unsigned reg_no = 0; reg_no < 16; reg_no++) {
            reg32_t src_val = randgen32->pick_val(0);
            reg32_t dst_val = randgen32->pick_val(0);
            GenRegArg src_reg = GenRegArg(cpu, mem, reg_no, src_val);
            SpecRegArg dst_reg = SpecRegArg(cpu, mem, dst_name,
                                            dstp, dst_val, src_val);
            failure = failure ||
                do_binary_reg_reg<GenRegArg, SpecRegArg, DefValidationFunc>(
                    cpu, mem, "LDS", src_reg, dst_reg, DefValidationFunc());
        }
        return failure;
    }

    // LDS Rm, MACH
    // 0100mmmm00001010
    static int binary_lds_gen_mach(Sh4 *cpu, Memory *mem,
                                   RandGen32 *randgen32) {
        return binary_lds_test_wrapper(cpu, mem, randgen32,
                                       "MACH", &cpu->reg.mach);
    }

    // LDS Rm, MACL
    // 0100mmmm00011010
    static int binary_lds_gen_macl(Sh4 *cpu, Memory *mem,
                                   RandGen32 *randgen32) {
        return binary_lds_test_wrapper(cpu, mem, randgen32,
                                       "MACL", &cpu->reg.macl);
    }

    // LDS Rm, PR
    // 0100mmmm00101010
    static int binary_lds_gen_pr(Sh4 *cpu, Memory *mem,
                                   RandGen32 *randgen32) {
        return binary_lds_test_wrapper(cpu, mem, randgen32,
                                       "PR", &cpu->reg.pr);
    }

    static int binary_sts_test_wrapper(Sh4 *cpu, Memory *mem,
                                       RandGen32 *randgen32,
                                       char const *src_name, reg32_t *srcp) {
        int failure = 0;
        for (unsigned reg_no = 0; reg_no < 16; reg_no++) {
            reg32_t src_val = randgen32->pick_val(0);
            reg32_t dst_val = randgen32->pick_val(0);
            SpecRegArg src_reg = SpecRegArg(cpu, mem, src_name,
                                            srcp, src_val);
            GenRegArg dst_reg = GenRegArg(cpu, mem, reg_no, dst_val, src_val);
            failure = failure ||
                do_binary_reg_reg<SpecRegArg, GenRegArg, DefValidationFunc>(
                    cpu, mem, "STS", src_reg, dst_reg, DefValidationFunc());
        }
        return failure;
    }

    // STS MACH, Rn
    // 0000nnnn00001010
    static int binary_sts_mach_gen(Sh4 *cpu, Memory *mem,
                                   RandGen32 *randgen32) {
        return binary_sts_test_wrapper(cpu, mem, randgen32,
                                       "MACH", &cpu->reg.mach);
    }

    // STS MACL, Rn
    // 0000nnnn00011010
    static int binary_sts_macl_gen(Sh4 *cpu, Memory *mem,
                                   RandGen32 *randgen32) {
        return binary_sts_test_wrapper(cpu, mem, randgen32,
                                       "MACL", &cpu->reg.macl);
    }

    // STS PR, Rn
    // 0000nnnn00101010
    static int binary_sts_pr_gen(Sh4 *cpu, Memory *mem,
                                   RandGen32 *randgen32) {
        return binary_sts_test_wrapper(cpu, mem, randgen32,
                                       "PR", &cpu->reg.pr);
    }

    // LDS.L @Rm+, MACH
    // 0100mmmm00000110
    static int do_binary_ldsl_indgeninc_mach(Sh4 *cpu, Memory *mem,
                                             unsigned reg_no, addr32_t addr,
                                             uint32_t val) {
        Sh4Prog test_prog;
        std::stringstream ss;
        std::string cmd;

        ss << "LDS.L @R" << reg_no << "+, MACH\n";
        cmd = ss.str();
        test_prog.assemble(cmd);
        const Sh4Prog::InstList& inst = test_prog.get_prog();
        mem->load_program(0, inst.begin(), inst.end());

        reset_cpu(cpu);

        *cpu->gen_reg(reg_no) = addr;
        cpu->write_mem(&val, addr, sizeof(val));

        cpu->exec_inst();

        if (cpu->reg.mach != val || *cpu->gen_reg(reg_no) != (addr + 4)) {
            std::cout << "ERROR while running " << cmd << std::endl;
            std::cout << "expected val is " << std::hex << val << std::endl;
            std::cout << "actual val is " << cpu->reg.mach << std::endl;
            std::cout << "input addr is " << addr << std::endl;
            std::cout << "output addr is " << (addr + 4) << std::endl;
            return 1;
        }

        return 0;
    }

    static int binary_ldsl_indgeninc_mach(Sh4 *cpu, Memory *mem,
                                          RandGen32 *randgen32) {
        int failure = 0;

        for (unsigned reg_no = 0; reg_no < 16; reg_no++) {
            addr32_t addr = randgen32->pick_range(0, mem->get_size() - 5);
            uint32_t val = randgen32->pick_val(0);
            failure = failure ||
                do_binary_ldsl_indgeninc_mach(cpu, mem, reg_no, addr, val);
        }

        return failure;
    }

    // LDS.L @Rm+, MACL
    // 0100mmmm00010110
    static int do_binary_ldsl_indgeninc_macl(Sh4 *cpu, Memory *mem,
                                             unsigned reg_no, addr32_t addr,
                                             uint32_t val) {
        Sh4Prog test_prog;
        std::stringstream ss;
        std::string cmd;

        ss << "LDS.L @R" << reg_no << "+, MACL\n";
        cmd = ss.str();
        test_prog.assemble(cmd);
        const Sh4Prog::InstList& inst = test_prog.get_prog();
        mem->load_program(0, inst.begin(), inst.end());

        reset_cpu(cpu);

        *cpu->gen_reg(reg_no) = addr;
        cpu->write_mem(&val, addr, sizeof(val));

        cpu->exec_inst();

        if (cpu->reg.macl != val || *cpu->gen_reg(reg_no) != (addr + 4)) {
            std::cout << "ERROR while running " << cmd << std::endl;
            std::cout << "expected val is " << std::hex << val << std::endl;
            std::cout << "actual val is " << cpu->reg.macl << std::endl;
            std::cout << "input addr is " << addr << std::endl;
            std::cout << "output addr is " << (addr + 4) << std::endl;
            return 1;
        }

        return 0;
    }

    static int binary_ldsl_indgeninc_macl(Sh4 *cpu, Memory *mem,
                                          RandGen32 *randgen32) {
        int failure = 0;

        for (unsigned reg_no = 0; reg_no < 16; reg_no++) {
            addr32_t addr = randgen32->pick_range(0, mem->get_size() - 5);
            uint32_t val = randgen32->pick_val(0);
            failure = failure ||
                do_binary_ldsl_indgeninc_macl(cpu, mem, reg_no, addr, val);
        }

        return failure;
    }

    // LDS.L @Rm+, PR
    // 0100mmmm00100110
    static int do_binary_ldsl_indgeninc_pr(Sh4 *cpu, Memory *mem,
                                           unsigned reg_no, addr32_t addr,
                                           uint32_t val) {
        Sh4Prog test_prog;
        std::stringstream ss;
        std::string cmd;

        ss << "LDS.L @R" << reg_no << "+, PR\n";
        cmd = ss.str();
        test_prog.assemble(cmd);
        const Sh4Prog::InstList& inst = test_prog.get_prog();
        mem->load_program(0, inst.begin(), inst.end());

        reset_cpu(cpu);

        *cpu->gen_reg(reg_no) = addr;
        cpu->write_mem(&val, addr, sizeof(val));

        cpu->exec_inst();

        if (cpu->reg.pr != val || *cpu->gen_reg(reg_no) != (addr + 4)) {
            std::cout << "ERROR while running " << cmd << std::endl;
            std::cout << "expected val is " << std::hex << val << std::endl;
            std::cout << "actual val is " << cpu->reg.pr << std::endl;
            std::cout << "input addr is " << addr << std::endl;
            std::cout << "output addr is " << (addr + 4) << std::endl;
            return 1;
        }

        return 0;
    }

    static int binary_ldsl_indgeninc_pr(Sh4 *cpu, Memory *mem,
                                        RandGen32 *randgen32) {
        int failure = 0;

        for (unsigned reg_no = 0; reg_no < 16; reg_no++) {
            addr32_t addr = randgen32->pick_range(0, mem->get_size() - 5);
            uint32_t val = randgen32->pick_val(0);
            failure = failure ||
                do_binary_ldsl_indgeninc_pr(cpu, mem, reg_no, addr, val);
        }

        return failure;
    }

    // STS.L MACH, @-Rn
    // 0100mmmm00000010
    static int do_binary_stsl_mach_inddecgen(Sh4 *cpu, Memory *mem,
                                             unsigned reg_no, reg32_t mach_val,
                                             addr32_t addr) {
        Sh4Prog test_prog;
        std::stringstream ss;
        std::string cmd;

        ss << "STS.L MACH, @-R" << reg_no << "\n";
        cmd = ss.str();
        test_prog.assemble(cmd);
        const Sh4Prog::InstList& inst = test_prog.get_prog();
        mem->load_program(0, inst.begin(), inst.end());

        reset_cpu(cpu);

        *cpu->gen_reg(reg_no) = addr;
        cpu->reg.mach = mach_val;
        cpu->exec_inst();

        uint32_t mem_val;
        cpu->read_mem(&mem_val, addr - 4, sizeof(mem_val));

        if (mem_val != mach_val || *cpu->gen_reg(reg_no) != (addr - 4)) {
            std::cout << "ERROR while running " << cmd << std::endl;
            std::cout << "expected val is " << std::hex << mach_val << std::endl;
            std::cout << "actual val is " << mem_val << std::endl;
            std::cout << "input addr is " << addr << std::endl;
            std::cout << "output addr is " << (addr - 4) << std::endl;
            return 1;
        }

        return 0;
    }

    static int binary_stsl_mach_inddecgen(Sh4 *cpu, Memory *mem,
                                          RandGen32 *randgen32) {
        int failure = 0;

        for (unsigned reg_no = 0; reg_no < 16; reg_no++) {
            addr32_t addr = randgen32->pick_range(4, mem->get_size() - 1);
            reg32_t mach_val = randgen32->pick_val(0);
            failure = failure ||
                do_binary_stsl_mach_inddecgen(cpu, mem, reg_no,
                                              mach_val, addr);
        }

        return failure;
    }

    // STS.L MACL, @-Rn
    // 0100mmmm00010010
    static int do_binary_stsl_macl_inddecgen(Sh4 *cpu, Memory *mem,
                                             unsigned reg_no, reg32_t macl_val,
                                             addr32_t addr) {
        Sh4Prog test_prog;
        std::stringstream ss;
        std::string cmd;

        ss << "STS.L MACL, @-R" << reg_no << "\n";
        cmd = ss.str();
        test_prog.assemble(cmd);
        const Sh4Prog::InstList& inst = test_prog.get_prog();
        mem->load_program(0, inst.begin(), inst.end());

        reset_cpu(cpu);

        *cpu->gen_reg(reg_no) = addr;
        cpu->reg.macl = macl_val;
        cpu->exec_inst();

        uint32_t mem_val;
        cpu->read_mem(&mem_val, addr - 4, sizeof(mem_val));

        if (mem_val != macl_val || *cpu->gen_reg(reg_no) != (addr - 4)) {
            std::cout << "ERROR while running " << cmd << std::endl;
            std::cout << "expected val is " << std::hex << macl_val << std::endl;
            std::cout << "actual val is " << mem_val << std::endl;
            std::cout << "input addr is " << addr << std::endl;
            std::cout << "output addr is " << (addr - 4) << std::endl;
            return 1;
        }

        return 0;
    }

    static int binary_stsl_macl_inddecgen(Sh4 *cpu, Memory *mem,
                                          RandGen32 *randgen32) {
        int failure = 0;

        for (unsigned reg_no = 0; reg_no < 16; reg_no++) {
            addr32_t addr = randgen32->pick_range(4, mem->get_size() - 1);
            reg32_t macl_val = randgen32->pick_val(0);
            failure = failure ||
                do_binary_stsl_macl_inddecgen(cpu, mem, reg_no,
                                              macl_val, addr);
        }

        return failure;
    }

    // STS.L PR, @-Rn
    // 0100nnnn00100010
    static int do_binary_stsl_pr_inddecgen(Sh4 *cpu, Memory *mem,
                                           unsigned reg_no, reg32_t pr_val,
                                           addr32_t addr) {
        Sh4Prog test_prog;
        std::stringstream ss;
        std::string cmd;

        ss << "STS.L PR, @-R" << reg_no << "\n";
        cmd = ss.str();
        test_prog.assemble(cmd);
        const Sh4Prog::InstList& inst = test_prog.get_prog();
        mem->load_program(0, inst.begin(), inst.end());

        reset_cpu(cpu);

        *cpu->gen_reg(reg_no) = addr;
        cpu->reg.pr = pr_val;
        cpu->exec_inst();

        uint32_t mem_val;
        cpu->read_mem(&mem_val, addr - 4, sizeof(mem_val));

        if (mem_val != pr_val || *cpu->gen_reg(reg_no) != (addr - 4)) {
            std::cout << "ERROR while running " << cmd << std::endl;
            std::cout << "expected val is " << std::hex << pr_val << std::endl;
            std::cout << "actual val is " << mem_val << std::endl;
            std::cout << "input addr is " << addr << std::endl;
            std::cout << "output addr is " << (addr - 4) << std::endl;
            return 1;
        }

        return 0;
    }

    static int binary_stsl_pr_inddecgen(Sh4 *cpu, Memory *mem,
                                        RandGen32 *randgen32) {
        int failure = 0;

        for (unsigned reg_no = 0; reg_no < 16; reg_no++) {
            addr32_t addr = randgen32->pick_range(4, mem->get_size() - 1);
            reg32_t pr_val = randgen32->pick_val(0);
            failure = failure ||
                do_binary_stsl_pr_inddecgen(cpu, mem, reg_no,
                                              pr_val, addr);
        }

        return failure;
    }

    // CMP/PZ Rn
    // 0100nnnn00010001
    static int do_unary_cmppz_gen(Sh4 *cpu, Memory *mem, unsigned reg_no,
                                  int32_t reg_val) {
        Sh4Prog test_prog;
        std::stringstream ss;
        std::string cmd;

        ss << "CMP/PZ R" << reg_no << "\n";
        cmd = ss.str();
        test_prog.assemble(cmd);
        const Sh4Prog::InstList& inst = test_prog.get_prog();
        mem->load_program(0, inst.begin(), inst.end());

        reset_cpu(cpu);

        *cpu->gen_reg(reg_no) = reg_val;
        cpu->exec_inst();

        bool t_expect = reg_val >= 0;

        if (t_expect && !(cpu->reg.sr & Sh4::SR_FLAG_T_MASK) ||
            !t_expect && (cpu->reg.sr & Sh4::SR_FLAG_T_MASK)) {
            std::cout << "ERROR while running " << cmd << std::endl;
            std::cout << "expected t val is " << t_expect << std::endl;
            std::cout << "actual val is " <<
                bool(cpu->reg.sr & Sh4::SR_FLAG_T_MASK) << std::endl;
            std::cout << "input val is " << reg_val << std::endl;
            return 1;
        }

        return 0;
    }

    static int unary_cmppz_gen(Sh4 *cpu, Memory *mem, RandGen32 *randgen32) {
        int failure = 0;

        for (unsigned reg_no = 0; reg_no < 16; reg_no++) {
            failure = failure ||
                do_unary_cmppz_gen(cpu, mem, reg_no, randgen32->pick_val(0));
        }

        return failure;
    }

    // CMP/PL Rn
    // 0100nnnn00010101
    static int do_unary_cmppl_gen(Sh4 *cpu, Memory *mem, unsigned reg_no,
                                  int32_t reg_val) {
        Sh4Prog test_prog;
        std::stringstream ss;
        std::string cmd;

        ss << "CMP/PL R" << reg_no << "\n";
        cmd = ss.str();
        test_prog.assemble(cmd);
        const Sh4Prog::InstList& inst = test_prog.get_prog();
        mem->load_program(0, inst.begin(), inst.end());

        reset_cpu(cpu);

        *cpu->gen_reg(reg_no) = reg_val;
        cpu->exec_inst();

        bool t_expect = reg_val > 0;

        if (t_expect && !(cpu->reg.sr & Sh4::SR_FLAG_T_MASK) ||
            !t_expect && (cpu->reg.sr & Sh4::SR_FLAG_T_MASK)) {
            std::cout << "ERROR while running " << cmd << std::endl;
            std::cout << "expected t val is " << t_expect << std::endl;
            std::cout << "actual val is " <<
                bool(cpu->reg.sr & Sh4::SR_FLAG_T_MASK) << std::endl;
            std::cout << "input val is " << reg_val << std::endl;
            return 1;
        }

        return 0;
    }

    static int unary_cmppl_gen(Sh4 *cpu, Memory *mem, RandGen32 *randgen32) {
        int failure = 0;

        for (unsigned reg_no = 0; reg_no < 16; reg_no++) {
            failure = failure ||
                do_unary_cmppl_gen(cpu, mem, reg_no, randgen32->pick_val(0));
        }

        return failure;
    }

    // CMP/EQ #imm, R0
    // 10001000iiiiiiii
    static int do_binary_cmpeq_imm_gen(Sh4 *cpu, Memory *mem,
                                       uint8_t imm_val, reg32_t r0_val) {
        GenRegArg dst_reg(cpu, mem, 0, r0_val);
        ImmedArg src_val(cpu, mem, imm_val);
        bool t_expect = (r0_val == int32_t(int8_t(imm_val)));
        return do_binary_reg_reg(cpu, mem, "CMP/EQ", src_val, dst_reg,
                                 ValidateFlagT(t_expect));
    }

    static int binary_cmpeq_imm_gen(Sh4 *cpu, Memory *mem,
                                    RandGen32 *randgen32) {
        int failure = 0;

        for (unsigned val = 0; val <= 255; val++) {
            uint8_t imm_val = val;
            // test equality
            failure = failure ||
                do_binary_cmpeq_imm_gen(cpu, mem, imm_val, imm_val);

            // test inequality
            failure = failure ||
                do_binary_cmpeq_imm_gen(cpu, mem, imm_val,
                                        int32_t(int8_t(imm_val)));
        }

        return failure;
    }

    // CMP/EQ Rm, Rn
    // 0011nnnnmmmm0000
    static int do_binary_cmpeq_gen_gen(Sh4 *cpu, Memory *mem,
                                       unsigned reg1, unsigned reg2,
                                       reg32_t reg1_val, reg32_t reg2_val) {
        GenRegArg dst_reg(cpu, mem, reg2, reg2_val);
        GenRegArg src_reg(cpu, mem, reg1, reg1_val);
        bool t_expect = (reg2_val == reg1_val);
        return do_binary_reg_reg(cpu, mem, "CMP/EQ", src_reg, dst_reg,
                                 ValidateFlagT(t_expect));
    }

    static int binary_cmpeq_gen_gen(Sh4 *cpu, Memory *mem,
                                    RandGen32 *randgen32) {
        int failure = 0;

        for (unsigned reg1 = 0; reg1 < 16; reg1++) {
            reg32_t val1 = randgen32->pick_val(0);
            for (unsigned reg2 = 0; reg2 < 16; reg2++) {
                reg32_t val2;
                if (reg1 == reg2)
                    val2 = val1;
                else
                    val2 = randgen32->pick_val(0);

                // test equality
                failure = failure ||
                    do_binary_cmpeq_gen_gen(cpu, mem, reg1, reg2, val2, val2);

                // test (probable) inequality
                failure = failure ||
                    do_binary_cmpeq_gen_gen(cpu, mem, reg1, reg2, val1, val2);
            }
        }

        return failure;
    }

    // CMP/HS Rm, Rn
    // 0011nnnnmmmm0010
    static int do_binary_cmphs_gen_gen(Sh4 *cpu, Memory *mem,
                                       unsigned reg1, unsigned reg2,
                                       reg32_t reg1_val, reg32_t reg2_val) {
        GenRegArg dst_reg(cpu, mem, reg2, reg2_val);
        GenRegArg src_reg(cpu, mem, reg1, reg1_val);
        bool t_expect = (reg2_val >= reg1_val);
        return do_binary_reg_reg(cpu, mem, "CMP/HS", src_reg, dst_reg,
                                 ValidateFlagT(t_expect));
    }

    static int binary_cmphs_gen_gen(Sh4 *cpu, Memory *mem,
                                    RandGen32 *randgen32) {
        int failure = 0;

        for (unsigned reg1 = 0; reg1 < 16; reg1++) {
            for (unsigned reg2 = 0; reg2 < 16; reg2++) {
                reg32_t val2 = randgen32->pick_val(0);

                // test equality
                failure = failure ||
                    do_binary_cmphs_gen_gen(cpu, mem, reg1, reg2, val2, val2);

                if (reg1 != reg2) {
                    // test val1 < val2
                    failure = failure ||
                        do_binary_cmphs_gen_gen(cpu, mem, reg1, reg2,
                                                val2 - 1, val2);

                    // test val1 > val2
                    failure = failure ||
                        do_binary_cmphs_gen_gen(cpu, mem, reg1, reg2,
                                                val2 + 1, val2);
                }
            }
        }

        return failure;
    }

    // CMP/GE Rm, Rn
    // 0011nnnnmmmm0011
    static int do_binary_cmpge_gen_gen(Sh4 *cpu, Memory *mem,
                                       unsigned reg1, unsigned reg2,
                                       reg32_t reg1_val, reg32_t reg2_val) {
        GenRegArg dst_reg(cpu, mem, reg2, reg2_val);
        GenRegArg src_reg(cpu, mem, reg1, reg1_val);
        bool t_expect = (int32_t(reg2_val) >= int32_t(reg1_val));
        return do_binary_reg_reg(cpu, mem, "CMP/GE", src_reg, dst_reg,
                                 ValidateFlagT(t_expect));
    }

    static int binary_cmpge_gen_gen(Sh4 *cpu, Memory *mem,
                                    RandGen32 *randgen32) {
        int failure = 0;

        for (unsigned reg1 = 0; reg1 < 16; reg1++) {
            for (unsigned reg2 = 0; reg2 < 16; reg2++) {
                reg32_t val2 = randgen32->pick_val(0);

                // test equality
                failure = failure ||
                    do_binary_cmpge_gen_gen(cpu, mem, reg1, reg2, val2, val2);

                if (reg1 != reg2) {
                    // test val1 < val2
                    failure = failure ||
                        do_binary_cmpge_gen_gen(cpu, mem, reg1, reg2,
                                                val2 - 1, val2);

                    // test val1 > val2
                    failure = failure ||
                        do_binary_cmpge_gen_gen(cpu, mem, reg1, reg2,
                                                val2 + 1, val2);
                }
            }
        }

        return failure;
    }

    // CMP/HI Rm, Rn
    // 0011nnnnmmmm0110
    static int do_binary_cmphi_gen_gen(Sh4 *cpu, Memory *mem,
                                       unsigned reg1, unsigned reg2,
                                       reg32_t reg1_val, reg32_t reg2_val) {
        GenRegArg dst_reg(cpu, mem, reg2, reg2_val);
        GenRegArg src_reg(cpu, mem, reg1, reg1_val);
        bool t_expect = (reg2_val > reg1_val);
        return do_binary_reg_reg(cpu, mem, "CMP/HI", src_reg, dst_reg,
                                 ValidateFlagT(t_expect));
    }

    static int binary_cmphi_gen_gen(Sh4 *cpu, Memory *mem,
                                    RandGen32 *randgen32) {
        int failure = 0;

        for (unsigned reg1 = 0; reg1 < 16; reg1++) {
            for (unsigned reg2 = 0; reg2 < 16; reg2++) {
                reg32_t val2 = randgen32->pick_val(0);

                // test equality
                failure = failure ||
                    do_binary_cmphi_gen_gen(cpu, mem, reg1, reg2, val2, val2);

                if (reg1 != reg2) {
                    // test val1 < val2
                    failure = failure ||
                        do_binary_cmphi_gen_gen(cpu, mem, reg1, reg2,
                                                val2 - 1, val2);

                    // test val1 > val2
                    failure = failure ||
                        do_binary_cmphi_gen_gen(cpu, mem, reg1, reg2,
                                                val2 + 1, val2);
                }
            }
        }

        return failure;
    }

    // CMP/GT Rm, Rn
    // 0011nnnnmmmm0111
    static int do_binary_cmpgt_gen_gen(Sh4 *cpu, Memory *mem,
                                       unsigned reg1, unsigned reg2,
                                       reg32_t reg1_val, reg32_t reg2_val) {
        GenRegArg dst_reg(cpu, mem, reg2, reg2_val);
        GenRegArg src_reg(cpu, mem, reg1, reg1_val);
        bool t_expect = (int32_t(reg2_val) > int32_t(reg1_val));
        return do_binary_reg_reg(cpu, mem, "CMP/GT", src_reg, dst_reg,
                                 ValidateFlagT(t_expect));
    }

    static int binary_cmpgt_gen_gen(Sh4 *cpu, Memory *mem,
                                    RandGen32 *randgen32) {
        int failure = 0;

        for (unsigned reg1 = 0; reg1 < 16; reg1++) {
            for (unsigned reg2 = 0; reg2 < 16; reg2++) {
                reg32_t val2 = randgen32->pick_val(0);

                // test equality
                failure = failure ||
                    do_binary_cmpgt_gen_gen(cpu, mem, reg1, reg2, val2, val2);

                if (reg1 != reg2) {
                    // test val1 < val2
                    failure = failure ||
                        do_binary_cmpgt_gen_gen(cpu, mem, reg1, reg2,
                                                val2 - 1, val2);

                    // test val1 > val2
                    failure = failure ||
                        do_binary_cmpgt_gen_gen(cpu, mem, reg1, reg2,
                                                val2 + 1, val2);
                }
            }
        }

        return failure;
    }

    static int do_binary_cmpstr_gen_gen(Sh4 *cpu, Memory *mem,
                                        unsigned reg1, unsigned reg2,
                                        reg32_t reg1_val, reg32_t reg2_val) {
        GenRegArg dst_reg(cpu, mem, reg2, reg2_val);
        GenRegArg src_reg(cpu, mem, reg1, reg1_val);
        bool t_expect = false;
        uint32_t mask = 0xff;
        for (int i = 0; i < 4; i++)
            if ((reg1_val & (0xff << (i * 8))) ==
                (reg2_val & (0xff << (i * 8))))
                t_expect = true;
        return do_binary_reg_reg(cpu, mem, "CMP/STR", src_reg, dst_reg,
                                 ValidateFlagT(t_expect));
    }

    static int binary_cmpstr_gen_gen(Sh4 *cpu, Memory *mem,
                                     RandGen32 *randgen32) {
        int failure = 0;

        for (unsigned reg1 = 0; reg1 < 16; reg1++) {
            for (unsigned reg2 = 0; reg2 < 16; reg2++) {
                reg32_t val2 = randgen32->pick_val(0);

                // test equality
                failure = failure ||
                    do_binary_cmpstr_gen_gen(cpu, mem, reg1, reg2, val2, val2);

                if (reg1 != reg2) {
                    // test partial equality
                    reg32_t val_tmp =
                        val2 ^ ~(0xff << (8 * randgen32->pick_range(0, 3)));
                    failure = failure ||
                        do_binary_cmpstr_gen_gen(cpu, mem, reg1, reg2,
                                                 val2, val_tmp);

                    // test (probable) inequality
                    failure = failure ||
                        do_binary_cmpstr_gen_gen(cpu, mem, reg1, reg2,
                                                 val2, randgen32->pick_val(0));
                }
            }
        }

        return failure;
    }

    // TST Rm, Rn
    // 0010nnnnmmmm1000
    static int binary_tst_gen_gen(Sh4 *cpu, Memory *mem, RandGen32 *randgen32) {
        int failure = 0;

        for (unsigned reg1_no = 0; reg1_no < 16; reg1_no++) {
            for (unsigned reg2_no = 0; reg2_no < 16; reg2_no++) {
                reg32_t reg1_val = randgen32->pick_val(0);
                reg32_t reg2_val = reg1_val;

                if (reg1_no != reg2_no)
                    reg2_val = randgen32->pick_val(0);

                GenRegArg dst_reg(cpu, mem, reg2_no, reg2_val);
                GenRegArg src_reg(cpu, mem, reg1_no, reg1_val);
                bool t_expect = !(reg1_val & reg2_val);
                failure = failure ||
                    do_binary_reg_reg(cpu, mem, "TST", src_reg, dst_reg,
                                      ValidateFlagT(t_expect));
            }
        }
        return failure;
    }

    // TAS.B @Rn
    // 0100nnnn00011011
    static int do_unary_tasb_indgen(Sh4 *cpu, Memory *mem, unsigned reg_no,
                                    addr32_t addr, uint8_t val) {
        Sh4Prog test_prog;
        std::stringstream ss;
        std::string cmd;

        ss << "TAS.B @R" << reg_no << "\n";
        cmd = ss.str();
        test_prog.assemble(cmd);
        const Sh4Prog::InstList& inst = test_prog.get_prog();
        mem->load_program(0, inst.begin(), inst.end());

        reset_cpu(cpu);

        *cpu->gen_reg(reg_no) = addr;
        cpu->write_mem(&val, addr, sizeof(val));
        cpu->exec_inst();
        bool t_expect = !val;

        if (t_expect && !(cpu->reg.sr & Sh4::SR_FLAG_T_MASK) ||
            !t_expect && (cpu->reg.sr & Sh4::SR_FLAG_T_MASK)) {
            std::cout << "ERROR while running " << cmd << std::endl;
            std::cout << "expected t val is " << t_expect << std::endl;
            std::cout << "actual t val is " <<
                bool(cpu->reg.sr & Sh4::SR_FLAG_T_MASK) << std::endl;
            std::cout << "val is " << std::hex << val << std::endl;
            std::cout << "addr is " << addr << std::endl;
            return 1;
        }

        return 0;
    }

    static int unary_tasb_indgen(Sh4 *cpu, Memory *mem, RandGen32 *randgen32) {
        int failure = 0;
        for (unsigned reg_no = 0; reg_no < 16; reg_no++) {
            addr32_t addr = randgen32->pick_range(0, mem->get_size() - 1);
            uint8_t val = randgen32->pick_val(0);

            // make sure 0 gets tested
            failure = failure ||
                do_unary_tasb_indgen(cpu, mem, reg_no, addr, 0);

            failure = failure ||
                do_unary_tasb_indgen(cpu, mem, reg_no, addr, val);
        }
        return failure;
    }

    // TST #imm, R0
    // 11001000iiiiiiii
    static int binary_tst_imm_r0(Sh4 *cpu, Memory *mem, RandGen32 *randgen32) {
        int failure = 0;

        for (unsigned imm_val = 0; imm_val < 256; imm_val++) {
            reg32_t reg1_val = randgen32->pick_val(0);

            ImmedArg src_val(cpu, mem, imm_val);
            GenRegArg dst_reg(cpu, mem, 0, reg1_val);
            bool t_expect = !(reg1_val & imm_val);
            failure = failure ||
                do_binary_reg_reg(cpu, mem, "TST", src_val, dst_reg,
                                  ValidateFlagT(t_expect));
        }
        return failure;
    }

    // TST.B #imm, @(R0, GBR)
    // 11001100iiiiiiii
    static int binary_tstb_imm_ind_r0_gbr(Sh4 *cpu, Memory *mem,
                                          RandGen32 *randgen32) {
        int failure = 0;

        for (unsigned imm_val = 0; imm_val < 256; imm_val++) {
            reg32_t gbr_val = randgen32->pick_range(0, mem->get_size() - 1) / 2;
            reg32_t r0_val = randgen32->pick_range(0, mem->get_size() - 1) / 2;
            uint8_t mem_val = randgen32->pick_val(0);
            bool t_expect = !(mem_val & uint8_t(imm_val));

            ImmedArg immed_arg(cpu, mem, imm_val);
            GenRegArg r0_arg(cpu, mem, 0, r0_val);
            SpecRegArg gbr_arg(cpu, mem, "GBR", &cpu->reg.gbr, gbr_val);

            BinIndArg<GenRegArg, SpecRegArg, 1> addr_arg(cpu, mem, r0_arg,
                                                         gbr_arg, mem_val,
                                                         mem_val);

            failure = failure ||
                do_binary_reg_reg(cpu, mem, "TST.B", immed_arg, addr_arg,
                                  ValidateFlagT(t_expect));
        }

        return failure;
    }

    // AND Rm, Rn
    // 0010nnnnmmmm1001
    static int binary_and_gen_gen(Sh4 *cpu, Memory *mem, RandGen32 *randgen32) {
        int failure = 0;
        for (int reg1_no = 0; reg1_no <= 15; reg1_no++) {
            for (int reg2_no = 0; reg2_no <= 15; reg2_no++) {
                reg32_t src_val = randgen32->pick_val(0);
                reg32_t dst_val = (reg1_no == reg2_no ? src_val :
                                   randgen32->pick_val(0));
                GenRegArg src(cpu, mem, reg1_no, src_val,
                              src_val == dst_val ? src_val & dst_val : src_val);
                GenRegArg dst(cpu, mem, reg2_no, dst_val, src_val & dst_val);
                failure = failure || do_binary_reg_reg(cpu, mem, "AND",
                                                       src, dst,
                                                       DefValidationFunc());
            }
        }
        return failure;
    }

    // AND #imm, R0
    // 11001001iiiiiiii
    static int binary_and_imm_r0(Sh4 *cpu, Memory *mem, RandGen32 *randgen32) {
        reg32_t initial_val = randgen32->pick_val(0);
        int failure = 0;
        for (unsigned imm_val = 0; imm_val <= 0xff; imm_val++) {
            uint32_t dst_val = randgen32->pick_val(0);
            ImmedArg src_op(cpu, mem, imm_val);
            GenRegArg dst_op(cpu, mem, 0, dst_val,
                             dst_val & uint32_t(imm_val));
            failure = failure ||
                do_binary_reg_reg(cpu, mem, "AND", src_op, dst_op,
                                  DefValidationFunc());
        }
        return failure;
    }

    // AND.B #imm, @(R0, GBR)
    // 11001101iiiiiiii
    static int binary_andb_imm_binind_r0_gbr(Sh4 *cpu, Memory *mem,
                                             RandGen32 *randgen32) {
        int failure = 0;

        for (unsigned imm_val = 0; imm_val < 256; imm_val++) {
            reg32_t gbr_val = randgen32->pick_range(0, mem->get_size() - 1) / 2;
            reg32_t r0_val = randgen32->pick_range(0, mem->get_size() - 1) / 2;
            uint8_t mem_val = randgen32->pick_val(0);
            uint8_t val_expect = mem_val & uint8_t(imm_val);

            ImmedArg immed_arg(cpu, mem, imm_val);
            GenRegArg r0_arg(cpu, mem, 0, r0_val);
            SpecRegArg gbr_arg(cpu, mem, "GBR", &cpu->reg.gbr, gbr_val);

            BinIndArg<GenRegArg, SpecRegArg, 1> addr_arg(cpu, mem, r0_arg,
                                                         gbr_arg, mem_val,
                                                         val_expect);

            failure = failure ||
                do_binary_reg_reg(cpu, mem, "AND.B", immed_arg, addr_arg,
                                  DefValidationFunc());
        }

        return failure;
    }

    // OR Rm, Rn
    // 0010nnnnmmmm1011
    static int binary_or_gen_gen(Sh4 *cpu, Memory *mem, RandGen32 *randgen32) {
        int failure = 0;
        for (int reg1_no = 0; reg1_no <= 15; reg1_no++) {
            for (int reg2_no = 0; reg2_no <= 15; reg2_no++) {
                reg32_t src_val = randgen32->pick_val(0);
                reg32_t dst_val = (reg1_no == reg2_no ? src_val :
                                   randgen32->pick_val(0));
                GenRegArg src(cpu, mem, reg1_no, src_val,
                              src_val == dst_val ? src_val | dst_val : src_val);
                GenRegArg dst(cpu, mem, reg2_no, dst_val, src_val | dst_val);
                failure = failure || do_binary_reg_reg(cpu, mem, "OR",
                                                       src, dst,
                                                       DefValidationFunc());
            }
        }
        return failure;
    }

    // OR #imm, R0
    // 11001011iiiiiiii
    static int binary_or_imm_r0(Sh4 *cpu, Memory *mem, RandGen32 *randgen32) {
        reg32_t initial_val = randgen32->pick_val(0);
        int failure = 0;
        for (unsigned imm_val = 0; imm_val <= 0xff; imm_val++) {
            uint32_t dst_val = randgen32->pick_val(0);
            ImmedArg src_op(cpu, mem, imm_val);
            GenRegArg dst_op(cpu, mem, 0, dst_val,
                             dst_val | uint32_t(imm_val));
            failure = failure ||
                do_binary_reg_reg(cpu, mem, "OR", src_op, dst_op,
                                  DefValidationFunc());
        }
        return failure;
    }

    // OR.B #imm, @(R0, GBR)
    // 11001111iiiiiiii
    static int binary_orb_imm_binind_r0_gbr(Sh4 *cpu, Memory *mem,
                                            RandGen32 *randgen32) {
        int failure = 0;

        for (unsigned imm_val = 0; imm_val < 256; imm_val++) {
            reg32_t gbr_val = randgen32->pick_range(0, mem->get_size() - 1) / 2;
            reg32_t r0_val = randgen32->pick_range(0, mem->get_size() - 1) / 2;
            uint8_t mem_val = randgen32->pick_val(0);
            uint8_t val_expect = mem_val | uint8_t(imm_val);

            ImmedArg immed_arg(cpu, mem, imm_val);
            GenRegArg r0_arg(cpu, mem, 0, r0_val);
            SpecRegArg gbr_arg(cpu, mem, "GBR", &cpu->reg.gbr, gbr_val);

            BinIndArg<GenRegArg, SpecRegArg, 1> addr_arg(cpu, mem, r0_arg,
                                                         gbr_arg, mem_val,
                                                         val_expect);

            failure = failure ||
                do_binary_reg_reg(cpu, mem, "OR.B", immed_arg, addr_arg,
                                  DefValidationFunc());
        }

        return failure;
    }

    // XOR Rm, Rn
    // 0010nnnnmmmm1010
    static int binary_xor_gen_gen(Sh4 *cpu, Memory *mem, RandGen32 *randgen32) {
        int failure = 0;
        for (int reg1_no = 0; reg1_no <= 15; reg1_no++) {
            for (int reg2_no = 0; reg2_no <= 15; reg2_no++) {
                reg32_t src_val = randgen32->pick_val(0);
                reg32_t dst_val = (reg1_no == reg2_no ? src_val :
                                   randgen32->pick_val(0));
                GenRegArg src(cpu, mem, reg1_no, src_val,
                              src_val == dst_val ? src_val ^ dst_val : src_val);
                GenRegArg dst(cpu, mem, reg2_no, dst_val, src_val ^ dst_val);
                failure = failure || do_binary_reg_reg(cpu, mem, "XOR",
                                                       src, dst,
                                                       DefValidationFunc());
            }
        }
        return failure;
    }

    // XOR #imm, R0
    // 11001010iiiiiiii
    static int binary_xor_imm_r0(Sh4 *cpu, Memory *mem, RandGen32 *randgen32) {
        reg32_t initial_val = randgen32->pick_val(0);
        int failure = 0;
        for (unsigned imm_val = 0; imm_val <= 0xff; imm_val++) {
            uint32_t dst_val = randgen32->pick_val(0);
            ImmedArg src_op(cpu, mem, imm_val);
            GenRegArg dst_op(cpu, mem, 0, dst_val,
                             dst_val ^ uint32_t(imm_val));
            failure = failure ||
                do_binary_reg_reg(cpu, mem, "XOR", src_op, dst_op,
                                  DefValidationFunc());
        }
        return failure;
    }

    // XOR.B #imm, @(R0, GBR)
    // 11001110iiiiiiii
    static int binary_xorb_imm_binind_r0_gbr(Sh4 *cpu, Memory *mem,
                                             RandGen32 *randgen32) {
        int failure = 0;

        for (unsigned imm_val = 0; imm_val < 256; imm_val++) {
            reg32_t gbr_val = randgen32->pick_range(0, mem->get_size() - 1) / 2;
            reg32_t r0_val = randgen32->pick_range(0, mem->get_size() - 1) / 2;
            uint8_t mem_val = randgen32->pick_val(0);
            uint8_t val_expect = mem_val ^ uint8_t(imm_val);

            ImmedArg immed_arg(cpu, mem, imm_val);
            GenRegArg r0_arg(cpu, mem, 0, r0_val);
            SpecRegArg gbr_arg(cpu, mem, "GBR", &cpu->reg.gbr, gbr_val);

            BinIndArg<GenRegArg, SpecRegArg, 1> addr_arg(cpu, mem, r0_arg,
                                                         gbr_arg, mem_val,
                                                         val_expect);

            failure = failure ||
                do_binary_reg_reg(cpu, mem, "XOR.B", immed_arg, addr_arg,
                                  DefValidationFunc());
        }

        return failure;
    }

    // NOT Rm, Rn
    // 0110nnnnmmmm0111
    static int binary_not_gen_gen(Sh4 *cpu, Memory *mem, RandGen32 *randgen32) {
        int failure = 0;
        for (int reg1_no = 0; reg1_no <= 15; reg1_no++) {
            for (int reg2_no = 0; reg2_no <= 15; reg2_no++) {
                reg32_t src_val = randgen32->pick_val(0);
                reg32_t dst_val = (reg1_no == reg2_no ? src_val :
                                   randgen32->pick_val(0));
                GenRegArg src(cpu, mem, reg1_no, src_val,
                              src_val == dst_val ? ~src_val : src_val);
                GenRegArg dst(cpu, mem, reg2_no, dst_val, ~src_val);
                failure = failure || do_binary_reg_reg(cpu, mem, "NOT",
                                                       src, dst,
                                                       DefValidationFunc());
            }
        }
        return failure;
    }
};

struct inst_test {
    char const *name;
    inst_test_func_t func;
} inst_tests[] = {
    { "nop_test", &Sh4InstTests::nop_test },
    { "add_immed_test", &Sh4InstTests::add_immed_test },
    { "add_gen_gen_test", &Sh4InstTests::add_gen_gen_test },
    { "addc_gen_gen_test", &Sh4InstTests::addc_gen_gen_test },
    { "addv_gen_gen_test", &Sh4InstTests::addv_gen_gen_test },
    { "sub_gen_gen_test", &Sh4InstTests::sub_gen_gen_test },
    { "subc_gen_gen_test", &Sh4InstTests::subc_gen_gen_test },
    { "subv_gen_gen_test", &Sh4InstTests::subv_gen_gen_test },
    { "movt_unary_gen_test", &Sh4InstTests::movt_unary_gen_test },
    { "mov_binary_imm_gen_test", &Sh4InstTests::mov_binary_imm_gen_test },
    { "movw_binary_binind_disp_pc_gen",
      &Sh4InstTests::movw_binary_binind_disp_pc_gen },
    { "movl_binary_binind_disp_pc_gen",
      &Sh4InstTests::movl_binary_binind_disp_pc_gen },
    { "mov_binary_gen_gen",
      &Sh4InstTests::mov_binary_gen_gen },
    { "movb_binary_gen_indgen", &Sh4InstTests::movb_binary_gen_indgen },
    { "movw_binary_gen_indgen", &Sh4InstTests::movw_binary_gen_indgen },
    { "movl_binary_gen_indgen", &Sh4InstTests::movl_binary_gen_indgen },
    { "movb_binary_indgen_gen", &Sh4InstTests::movb_binary_indgen_gen },
    { "movw_binary_indgen_gen", &Sh4InstTests::movw_binary_indgen_gen },
    { "movl_binary_indgen_gen", &Sh4InstTests::movl_binary_indgen_gen },
    { "movb_binary_gen_inddecgen", &Sh4InstTests::movb_binary_gen_inddecgen },
    { "movw_binary_gen_inddecgen", &Sh4InstTests::movw_binary_gen_inddecgen },
    { "movl_binary_gen_inddecgen", &Sh4InstTests::movl_binary_gen_inddecgen },
    { "movb_binary_indgeninc_gen", &Sh4InstTests::movb_binary_indgeninc_gen },
    { "movw_binary_indgeninc_gen", &Sh4InstTests::movw_binary_indgeninc_gen },
    { "movl_binary_indgeninc_gen", &Sh4InstTests::movl_binary_indgeninc_gen },
    { "movb_binary_r0_binind_disp_gen",
      &Sh4InstTests::movb_binary_r0_binind_disp_gen },
    { "movw_binary_r0_binind_disp_gen",
      &Sh4InstTests::movw_binary_r0_binind_disp_gen },
    { "movl_binary_gen_binind_disp_gen",
      &Sh4InstTests::movl_binary_gen_binind_disp_gen },
    { "movb_binary_binind_disp_gen_r0",
      &Sh4InstTests::movb_binary_binind_disp_gen_r0 },
    { "movw_binary_binind_disp_gen_r0",
      &Sh4InstTests::movw_binary_binind_disp_gen_r0 },
    { "movl_binary_binind_disp_gen_gen",
      &Sh4InstTests::movl_binary_binind_disp_gen_gen },
    { "movb_gen_binind_r0_gen", &Sh4InstTests::movb_gen_binind_r0_gen },
    { "movw_gen_binind_r0_gen", &Sh4InstTests::movw_gen_binind_r0_gen },
    { "movl_gen_binind_r0_gen", &Sh4InstTests::movl_gen_binind_r0_gen },
    { "binary_movb_binind_r0_gen_gen",
      &Sh4InstTests::binary_movb_binind_r0_gen_gen },
    { "binary_movw_binind_r0_gen_gen",
      &Sh4InstTests::binary_movw_binind_r0_gen_gen },
    { "binary_movl_binind_r0_gen_gen",
      &Sh4InstTests::binary_movl_binind_r0_gen_gen },
    { "binary_movb_r0_binind_disp_gbr",
      &Sh4InstTests::binary_movb_r0_binind_disp_gbr },
    { "binary_movw_r0_binind_disp_gbr",
      &Sh4InstTests::binary_movw_r0_binind_disp_gbr },
    { "binary_movl_r0_binind_disp_gbr",
      &Sh4InstTests::binary_movl_r0_binind_disp_gbr },
    { "binary_movb_binind_disp_gbr_r0",
      &Sh4InstTests::binary_movb_binind_disp_gbr_r0 },
    { "binary_movw_binind_disp_gbr_r0",
      &Sh4InstTests::binary_movw_binind_disp_gbr_r0 },
    { "binary_movl_binind_disp_gbr_r0",
      &Sh4InstTests::binary_movl_binind_disp_gbr_r0 },
    { "binary_mova_binind_disp_pc_r0",
      &Sh4InstTests::binary_mova_binind_disp_pc_r0 },
    { "binary_ldc_gen_sr", &Sh4InstTests::binary_ldc_gen_sr },
    { "binary_ldc_gen_gbr", &Sh4InstTests::binary_ldc_gen_gbr },
    { "binary_ldc_gen_vbr", &Sh4InstTests::binary_ldc_gen_vbr },
    { "binary_ldc_gen_ssr", &Sh4InstTests::binary_ldc_gen_ssr },
    { "binary_ldc_gen_spc", &Sh4InstTests::binary_ldc_gen_spc },
    { "binary_ldc_gen_bank", &Sh4InstTests::binary_ldc_gen_bank },
    { "binary_ldcl_indgeninc_sr", &Sh4InstTests::binary_ldcl_indgeninc_sr },
    { "binary_ldcl_indgeninc_gbr", &Sh4InstTests::binary_ldcl_indgeninc_gbr },
    { "binary_ldcl_indgeninc_vbr", &Sh4InstTests::binary_ldcl_indgeninc_vbr },
    { "binary_ldcl_indgeninc_ssr", &Sh4InstTests::binary_ldcl_indgeninc_ssr },
    { "binary_ldcl_indgeninc_spc", &Sh4InstTests::binary_ldcl_indgeninc_spc },
    { "binary_ldcl_indgeninc_dbr", &Sh4InstTests::binary_ldcl_indgeninc_dbr },
    { "binary_stc_sr_gen", &Sh4InstTests::binary_stc_sr_gen },
    { "binary_stc_gbr_gen", &Sh4InstTests::binary_stc_gbr_gen },
    { "binary_stc_vbr_gen", &Sh4InstTests::binary_stc_vbr_gen },
    { "binary_stc_ssr_gen", &Sh4InstTests::binary_stc_ssr_gen },
    { "binary_stc_spc_gen", &Sh4InstTests::binary_stc_spc_gen },
    { "binary_stc_sgr_gen", &Sh4InstTests::binary_stc_sgr_gen },
    { "binary_stc_dbr_gen", &Sh4InstTests::binary_stc_dbr_gen },
    { "binary_stcl_sr_inddecgen", &Sh4InstTests::binary_stcl_sr_inddecgen },
    { "binary_stcl_gbr_inddecgen", &Sh4InstTests::binary_stcl_gbr_inddecgen },
    { "binary_stcl_vbr_inddecgen", &Sh4InstTests::binary_stcl_vbr_inddecgen },
    { "binary_stcl_ssr_inddecgen", &Sh4InstTests::binary_stcl_ssr_inddecgen },
    { "binary_stcl_spc_inddecgen", &Sh4InstTests::binary_stcl_spc_inddecgen },
    { "binary_stcl_sgr_inddecgen", &Sh4InstTests::binary_stcl_sgr_inddecgen },
    { "binary_stcl_dbr_inddecgen", &Sh4InstTests::binary_stcl_dbr_inddecgen },
    { "binary_ldcl_indgeninc_bank", &Sh4InstTests::binary_ldcl_indgeninc_bank },
    { "binary_stc_bank_gen", &Sh4InstTests::binary_stc_bank_gen },
    { "binary_stcl_bank_inddecgen", &Sh4InstTests::binary_stcl_bank_inddecgen },
    { "binary_lds_gen_mach", &Sh4InstTests::binary_lds_gen_mach },
    { "binary_lds_gen_macl", &Sh4InstTests::binary_lds_gen_macl },
    { "binary_lds_gen_pr", &Sh4InstTests::binary_lds_gen_pr },
    { "binary_sts_mach_gen", &Sh4InstTests::binary_sts_mach_gen },
    { "binary_sts_macl_gen", &Sh4InstTests::binary_sts_macl_gen },
    { "binary_sts_pr_gen", &Sh4InstTests::binary_sts_pr_gen },
    { "binary_ldsl_indgeninc_mach", &Sh4InstTests::binary_ldsl_indgeninc_mach },
    { "binary_ldsl_indgeninc_macl", &Sh4InstTests::binary_ldsl_indgeninc_macl },
    { "binary_ldsl_indgeninc_pr",  &Sh4InstTests::binary_ldsl_indgeninc_pr },
    { "binary_stsl_mach_inddecgen", &Sh4InstTests::binary_stsl_mach_inddecgen },
    { "binary_stsl_macl_inddecgen", &Sh4InstTests::binary_stsl_macl_inddecgen },
    { "binary_stsl_pr_inddecgen", &Sh4InstTests::binary_stsl_pr_inddecgen },
    { "unary_cmppz_gen", &Sh4InstTests::unary_cmppz_gen },
    { "unary_cmppl_gen", &Sh4InstTests::unary_cmppl_gen },
    { "binary_cmpeq_imm_gen", &Sh4InstTests::binary_cmpeq_imm_gen },
    { "binary_cmpeq_gen_gen", &Sh4InstTests::binary_cmpeq_gen_gen },
    { "binary_cmphs_gen_gen", &Sh4InstTests::binary_cmphs_gen_gen },
    { "binary_cmpge_gen_gen", &Sh4InstTests::binary_cmpge_gen_gen },
    { "binary_cmphi_gen_gen", &Sh4InstTests::binary_cmphi_gen_gen },
    { "binary_cmpgt_gen_gen", &Sh4InstTests::binary_cmpgt_gen_gen },
    { "binary_cmpstr_gen_gen", &Sh4InstTests::binary_cmpstr_gen_gen },
    { "binary_tst_gen_gen", &Sh4InstTests::binary_tst_gen_gen },
    { "unary_tasb_indgen", &Sh4InstTests::unary_tasb_indgen },
    { "binary_tst_imm_r0", &Sh4InstTests::binary_tst_imm_r0 },
    { "binary_tstb_imm_ind_r0_gbr", &Sh4InstTests::binary_tstb_imm_ind_r0_gbr },
    { "binary_and_gen_gen", &Sh4InstTests::binary_and_gen_gen },
    { "binary_and_imm_r0", &Sh4InstTests::binary_and_imm_r0 },
    { "binary_andb_imm_binind_r0_gbr", &Sh4InstTests::binary_andb_imm_binind_r0_gbr },
    { "binary_or_gen_gen", &Sh4InstTests::binary_or_gen_gen },
    { "binary_or_imm_r0", &Sh4InstTests::binary_or_imm_r0 },
    { "binary_orb_imm_binind_r0_gbr", &Sh4InstTests::binary_orb_imm_binind_r0_gbr },
    { "binary_xor_gen_gen", &Sh4InstTests::binary_xor_gen_gen },
    { "binary_xor_imm_r0", &Sh4InstTests::binary_xor_imm_r0 },
    { "binary_xorb_imm_binind_r0_gbr", &Sh4InstTests::binary_xorb_imm_binind_r0_gbr },
    { "binary_not_gen_gen", &Sh4InstTests::binary_not_gen_gen },
    { NULL }
};

int main(int argc, char **argv) {
    Memory mem(16 * 1024 * 1024);
    Sh4 cpu(&mem);
    struct inst_test *test = inst_tests;
    int n_success = 0, n_tests = 0;
    unsigned int seed = time(NULL);
    int opt;

    while ((opt = getopt(argc, argv, "s:")) > 0) {
        if (opt == 's')
            seed = atoi(optarg);
    }

    try {
        RandGen32 randgen32(seed);
        randgen32.reset();

        while (test->name) {
            std::cout << "Trying " << test->name << "..." << std::endl;

            int test_ret = test->func(&cpu, &mem, &randgen32);

            if (test_ret != 0)
                std::cout << test->name << " FAIL" << std::endl;
            else {
                std::cout << test->name << " SUCCESS" << std::endl;
                n_success++;
            }

            test++;
            n_tests++;
        }
    } catch (UnimplementedError excp) {
        std::cerr << "ERROR: " << excp.what() << std::endl;
        return 1;
    }

    double percent = 100.0 * double(n_success) / double(n_tests);
    std::cout << std::dec << n_tests << " tests run - " << n_success <<
        " successes " << "(" << percent << "%)" << std::endl;

    if (n_success == n_tests)
        return 0;
    return 1;
}
