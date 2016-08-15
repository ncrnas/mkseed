#include "gtest/gtest.h"
#include "mkseed_core.hpp"
#include "get_data_path.hpp"
//
//
//class fasta_test : public ::testing::Test
//{
//protected:
//    virtual void SetUp() {
//        std::string tmpf = std::tmpnam(nullptr);
//
//
//    }
//
//    virtual void TearDown() {
//        fclose(tmpf);
//
//    }
//private:
//    std::FILE* tmpf;
//};
//
//TEST_F(fasta_test, CaseDataWipe) {
//
//}

TEST(fast_check, mirna_fast) {
    seqan::CharString dfile = STRINGIZE(TEST_DATA_PATH);
    dfile += "test_mature_1.fa";

    std::cout << dfile << std::endl;

    std::fstream mirna_fa(seqan::toCString(dfile), std::ios::binary | std::ios::in);
    EXPECT_TRUE(mirna_fa.good());
}

TEST(fast_check, mrna_fast) {
    seqan::CharString dfile = STRINGIZE(TEST_DATA_PATH);
    dfile += "test_refseq_1.fa";

    std::fstream mrna_fa(seqan::toCString(dfile), std::ios::binary | std::ios::in);
    EXPECT_TRUE(mrna_fa.good());
}