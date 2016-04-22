require_relative "../read"

describe Read do
  describe "::k_plus_one_mers" do
    it "yields [kmer, start_pos]" do
      read = "abcde"
      ksize = 2

      expect { |b| Read.k_plus_one_mers(read, ksize, &b) }.
        to yield_successive_args(["abc", 0], ["bcd", 1], ["cde", 2])

    end
  end
end
