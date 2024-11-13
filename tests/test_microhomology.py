import unittest
from destruct.utils.misc import homology_consistent_breakpoint
import destruct.utils.misc


def create_sequence(chromosome, strand, position, reference_sequences):
    breakend_sequences = ['', '']
    expected_strands = ('+', '-')
    length = 10
    for side in (0, 1):
        if strand[side] == '+':
            start = position[side] - length + 1
            end = position[side]
        else:
            start = position[side]
            end = position[side] + length - 1
        breakend_sequences[side] = reference_sequences[chromosome[side]][start-1:end]
        if strand[side] != expected_strands[side]:
            breakend_sequences[side] = destruct.utils.misc.reverse_complement(breakend_sequences[side])
    return breakend_sequences[0] + breakend_sequences[1]


class TestMicrohomology(unittest.TestCase):
    
    def test_case_1(self):
        mock_genome = {   
            # ------------AGCTAGCTTACATGCTAGCT
            # ------------||||||||||
            "chr1": "AGCTTAGCTAGCTTAAGATCGATCGGGATCGTAGCTAGCTA",
            # -----------AGCTAGCTTACATGCTAGCT
            # ---------------------||||||||||
            "chr2": "CGTAGCTAGGCTAGCATGCTAGCTAGATCGTAGCTAACGTAG",
            "chr3": "TACGATCGTAGCATCGATCGGCTAGCTACGATCGTAGCATCG"
        }

        chromosome = ["chr1", "chr2"]
        strand = ["+", "-"]
        position = [15, 15]
        maxOffset = 10

        result = homology_consistent_breakpoint(chromosome, strand, position, mock_genome, maxOffset)
        self.assertEqual(result, 0)


    def test_case_2(self):
        mock_genome = {
            # ------------AGCTAGCTTACATGCTAGCT
            # ------------||||||||||X
            "chr1": "AGCTTAGCTAGCTTACGATCGATCGGGATCGTAGCTAGCTA",
            # -----------AGCTAGCTTACATGCTAGCT
            # ---------------------||||||||||
            "chr2": "CGTAGCTAGGCTAGCATGCTAGCTAGATCGTAGCTAACGTAG",
            "chr3": "TACGATCGTAGCATCGATCGGCTAGCTACGATCGTAGCATCG"
        }

        chromosome = ["chr1", "chr2"]
        strand = ["+", "-"]
        position = [15, 15]
        maxOffset = 10

        result = homology_consistent_breakpoint(chromosome, strand, position, mock_genome, maxOffset)
        self.assertEqual(result, 1)


    def test_case_3(self):
        mock_genome = {
            # ------------AGCTAGCTTACATGCTAGCT
            # ------------||||||||||X
            "chr1": "AGCTTAGCTAGCTTACGATCGATCGGGATCGTAGCTAGCTA",
            # -----------AGCTAGCTTACATGCTAGCT
            # --------------------X||||||||||
            "chr2": "CGTAGCTAGGCTAACATGCTAGCTAGATCGTAGCTAACGTAG",
            "chr3": "TACGATCGTAGCATCGATCGGCTAGCTACGATCGTAGCATCG"
        }

        chromosome = ["chr1", "chr2"]
        strand = ["+", "-"]
        position = [15, 15]
        maxOffset = 10

        result = homology_consistent_breakpoint(chromosome, strand, position, mock_genome, maxOffset)
        self.assertEqual(result, 2)


    def test_case_4(self):
        mock_genome = {
            # ------------AGCTAGCTTACATGCTAGCT
            # ------------||||||||||XXX
            "chr1": "AGCTTAGCTAGCTTACATTCGATCGGGATCGTAGCTAGCTA",
            # -----------AGCTAGCTTACATGCTAGCT
            # ----------------XXXXX||||||||||
            "chr2": "CGTAGCTAGGCTTACATGCTAGCTAGATCGTAGCTAACGTAG",
            "chr3": "TACGATCGTAGCATCGATCGGCTAGCTACGATCGTAGCATCG"
        }

        chromosome = ["chr1", "chr2"]
        strand = ["+", "-"]
        position = [15, 15]
        maxOffset = 10

        result = homology_consistent_breakpoint(chromosome, strand, position, mock_genome, maxOffset)
        self.assertEqual(result, 8)


    def test_case_5(self):
        mock_genome = {
            # ------------AGCTAGCTTAGTAAGCCTAG
            # ------------||||||||||
            "chr1": "AGCTTAGCTAGCTTACATTCGATCGGGATCGTAGCTAGCTA",
            # ------------CTAGGCTTACTAAGCTAGCT (rc)
            # ------------||||||||||
            "chr2": "CGTAGCTAGGCTTACATGCTAGCTAGATCGTAGCTAACGTAG",
            "chr3": "TACGATCGTAGCATCGATCGGCTAGCTACGATCGTAGCATCG"
        }

        chromosome = ["chr1", "chr2"]
        strand = ["+", "+"]
        position = [15, 15]
        maxOffset = 10

        result = homology_consistent_breakpoint(chromosome, strand, position, mock_genome, maxOffset)
        self.assertEqual(result, 0)


    def test_case_6(self):
        mock_genome = {
            # ------------AGCTAGCTTAGTAAGCCTAG
            # ------------||||||||||XXX
            "chr1": "AGCTTAGCTAGCTTAGTATTCGATCGGGATCGTAGCTAGCTA",
            # ------------CTAGGCTTACTAAGCTAGCT (rc)
            # ------------||||||||||
            "chr2": "CGTAGCTAGGCTTACATGCTAGCTAGATCGTAGCTAACGTAG",
            "chr3": "TACGATCGTAGCATCGATCGGCTAGCTACGATCGTAGCATCG"
        }

        chromosome = ["chr1", "chr2"]
        strand = ["+", "+"]
        position = [15, 15]
        maxOffset = 10

        result = homology_consistent_breakpoint(chromosome, strand, position, mock_genome, maxOffset)
        self.assertEqual(result, 3)


    def test_case_7(self):
        mock_genome = {
            # ------------AGCTAGCTTAGTAAGCCTAG
            # ------------||||||||||XXX
            "chr1": "AGCTTAGCTAGCTTAGTATTCGATCGGGATCGTAGCTAGCTA",
            # ------------CTAGGCTTACTAAGCTAGCT (rc)
            # ------------||||||||||XXX
            "chr2": "CGTAGCTAGGCTTACTAACTAGCTAGATCGTAGCTAACGTAG",
            "chr3": "TACGATCGTAGCATCGATCGGCTAGCTACGATCGTAGCATCG"
        }

        chromosome = ["chr1", "chr2"]
        strand = ["+", "+"]
        position = [15, 15]
        maxOffset = 10

        result = homology_consistent_breakpoint(chromosome, strand, position, mock_genome, maxOffset)
        self.assertEqual(result, 6)


    def test_case_8(self):
        mock_genome = {
            # -----------AGCTAGTTAGAGTATTCGAT
            # ---------------------||||||||||
            "chr1": "AGCTTAGCTAGCTTAGTATTCGATCGGGATCGTAGCTAGCTA",
            # -----------ATCGAATACTCTAACTAGCT
            # ---------------------||||||||||
            "chr2": "CGTAGCTAGGCTTACTAACTAGCTAGATCGTAGCTAACGTAG",
            "chr3": "TACGATCGTAGCATCGATCGGCTAGCTACGATCGTAGCATCG"
        }

        chromosome = ["chr1", "chr2"]
        strand = ["-", "-"]
        position = [15, 15]
        maxOffset = 10

        result = homology_consistent_breakpoint(chromosome, strand, position, mock_genome, maxOffset)
        self.assertEqual(result, 0)


    def test_case_9(self):
        mock_genome = {
            # -----------AGCTAGTTAGAGTATTCGAT
            # --------------------X||||||||||
            "chr1": "AGCTTAGCTAGCTGAGTATTCGATCGGGATCGTAGCTAGCTA",
            # -----------ATCGAATACTCTAACTAGCT
            # --------------------X||||||||||
            "chr2": "CGTAGCTAGGCTTTCTAACTAGCTAGATCGTAGCTAACGTAG",
            "chr3": "TACGATCGTAGCATCGATCGGCTAGCTACGATCGTAGCATCG"
        }

        chromosome = ["chr1", "chr2"]
        strand = ["-", "-"]
        position = [15, 15]
        maxOffset = 10
        
        result = homology_consistent_breakpoint(chromosome, strand, position, mock_genome, maxOffset)
        self.assertEqual(result, 2)


if __name__ == "__main__":
    unittest.main()
