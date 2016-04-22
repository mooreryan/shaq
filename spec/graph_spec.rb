# Copyright 2016 Ryan Moore
# Contact: moorer@udel.edu
#
# This file is part of Shaq.
#
# Shaq is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Shaq is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
# License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Shaq.  If not, see <http://www.gnu.org/licenses/>.

require_relative "../graph"

describe Graph do
  let(:ksize) { 3 }
  let(:graph) { Graph.new ksize }

  let(:header) { "read 1" }
  let(:read) { "aaataaac" }

  let(:db_graph) {
    { "aaa" => Set.new(["aat", "aac"]),
      "aat" => Set.new(["ata"]),
      "ata" => Set.new(["taa"]),
      "taa" => Set.new(["aaa"])}
  }

  let(:kmer_to_read) {
    { "aaa" => [["read 1", 0], ["read 1", 4]],
      "aat" => [["read 1", 1]],
      "ata" => [["read 1", 2]],
      "taa" => [["read 1", 3]],
      "aac" => [["read 1", 5]], }
  }

  let(:kmer_abundance) {
    { "aaa" => 2,
      "aat" => 1,
      "ata" => 1,
      "taa" => 1,
      "aac" => 1,}
  }

  let(:incounts) {
    { "aaa" => 1,
      "aat" => 1,
      "ata" => 1,
      "taa" => 1,
      "aac" => 1,}
  }


  before(:each) do
    graph.add_to_graph header, read
  end

  describe "#add_to_graph" do
    it "adds the read to the graph" do
      expect(graph.graph).to eq db_graph
    end

    it "adds the header and pos to the kmer_to_read hash" do
      expect(graph.kmer_to_read).to eq kmer_to_read
    end

    it "tracks number of incoming connections to kmers" do
      expect(graph.incounts).to eq incounts
    end
  end

  describe "#kmer_abundance" do
    it "returns a coll with kmer abundance" do
      expect(graph.kmer_abundance).to eq kmer_abundance
    end
  end
end
