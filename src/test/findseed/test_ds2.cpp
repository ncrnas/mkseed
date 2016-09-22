#include <iostream>
#include "test_ds.hpp"

namespace {

    class TestDS2 : public TestDS
    {
    protected:
        TestDS2() {
            IFNAME1 = (char *)"hsa_mirbase_1.txt";
            IFNAME2 = (char *)"utr2.txt";
            OFNAME1 = (char *)"test_output_2.txt";
            OFNAME2 = (char *)"test_ds2.txt";
        }
    };

    TEST_F(TestDS2, mrna_fasta) {
        read_files();

        EXPECT_EQ(1u, length(mrna_ids));

        const char *id1 = "hg18_refGene NM_005498 range=chr19:10544348-10544741 5'pad=0 3'pad=0 "
                "revComp=TRUE strand=- repeatMasking=none";
        EXPECT_STREQ(id1, seqan::toCString(mrna_ids[0]));

        EXPECT_EQ(1u, length(mrna_seqs));
        const char *seq1 = "AAGGGAGAAGAGAUGGGGGCUUGAACACGGGGCUUCCUUACAGCCCCGGA"
                "UGCAGAUUUUAGAGGGAGGGCAGGUGCGGGCUGUGUGUGUCUGUGUGAGG"
                "GCAGGUCCUGGACUUGGCAGUUUCUUGCUCCCAGCACCCGCCCCUUCCUC"
                "ACCUCUUCCUUAUUCCAUAGGCUGGGAGAGAAACUCUCUCUGCUUCCCUC"
                "GCCCUUGGAGCUUUCCCCAUCCCCCUGAUUUUAUAUGAAGAAAUAGAAGA"
                "GGGGCUUGAAGUCCUCCUCGCGAGUGCCUUCUUGCAAUUACCUGCCUUAG"
                "CGGGUGUUGCGGGUCCCUCCUUCACAGCCGCUGAGCCCAGAGGUCCCGCU"
                "GGCCCCUCCUCUGAAUUUUAGGAUGUCAUUAAAAAGAUGAAUCU";
        EXPECT_STREQ(seq1, seqan::toCString((seqan::CharString)mrna_seqs[0]));
    }

    TEST_F(TestDS2, comp_ngram) {
        comp_ngram();
    }

    TEST_F(TestDS2, comp_bruteforce) {
        comp_bruteforce();
    }

    TEST_F(TestDS2, comp_horspool) {
        comp_horspool();
    }

    TEST_F(TestDS2, comp_suffixarray) {
        comp_suffixarray();
    }
}