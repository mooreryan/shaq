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

require "parse_fasta"
require "set"

def branch_point? set
  set.count > 1
end

def debruijn kmer_info, in_out, revhash, hash, rec, k
  readname, read = rec
  # start_pos is zero based
  read.split("").each_cons(k+1).each_with_index do |kmer, start_pos|
    kmer = kmer.join

    k1 = kmer[0..-2]
    k2 = kmer[1..-1]

    k1_info = [readname, start_pos]
    # k2_info = [readname, start_pos+1]

    if kmer_info.has_key? k1
      kmer_info[k1] << k1_info
    else
      kmer_info[k1] = [k1_info]
    end

    # if kmer_info.has_key? k2
    #   kmer_info[k2] << k2_info
    # else
    #   kmer_info[k2] = [k2_info]
    # end

    if in_out.has_key? k1
      in_out[k1][:out] += 1
    else
      in_out[k1] = {in: 0, out: 1}
    end

    if in_out.has_key? k2
      in_out[k2][:in] += 1
    else
      in_out[k2] = {in: 1, out: 0}
    end

    if hash.has_key? k1
      hash[k1] << k2
    else
      hash[k1] = Set.new [k2]
    end

    if revhash.has_key? k2
      revhash[k2] << k1
    else
      revhash[k2] = Set.new [k1]
    end
  end

  [kmer_info, in_out, revhash, hash]
end

ksize = 9
hash = {}
revhash = {}
in_out = {}
kmer_info = {}

SeqFile.open(ARGV.first).each_record_fast do |head, seq|
  kmer_info, in_out, revhash, hash = debruijn kmer_info, in_out, revhash, hash, [head, seq], ksize
end

tmp_walk = ""
fin_walk = ""

contigs = {}
branch_point_kmers = {}

kmer = in_out.select { |kmer, connections| connections[:in].zero? }.first.first
in_out.delete kmer
tmp_walk = String.new kmer
next_kmer = ""
num = 0
while !next_kmer.nil?
  if hash.has_key?(kmer) && !hash[kmer].empty?
    if branch_point? hash[kmer]
      # this kmer will be visited again, mark it as being connected to
      # current contig
      if branch_point_kmers.has_key? kmer
        branch_point_kmers[kmer][:end_of] << num
        branch_point_kmers[kmer][:start_of] << num+1
      else
        branch_point_kmers[kmer] = { end_of: [num],
                                     start_of: [num+1] }
      end

      next_kmer = hash[kmer].first

      # check if there is a better kmer
      # this_readname, this_start_pos = kmer_info[kmer]

      # TODO this big giant loopy thing doesn't actually work cos all
      #   side by side kmers will be on the same read -- needs to go
      #   down a path and if it diverges from the read, go back and
      #   check other branch direction. something like that
      # break_loop = false
      # hash[kmer].each do |next_one|
      #   kmer_info[next_one].each do |next_readname, next_start_pos|
      #     kmer_info[kmer].each do |this_readname, this_start_pos|
      #       if this_readname == next_readname && this_start_pos == next_start_pos - 1
      #         next_kmer = next_one
      #         break_loop = true
      #         break
      #       end
      #     end

      #     break if break_loop
      #   end

      #   break if break_loop
      # end
      contigs[num] = tmp_walk
      num += 1
      tmp_walk = String.new( kmer[0] + next_kmer) # next_kmer # << next_kmer[-1]
    else
      next_kmer = hash[kmer].first
      tmp_walk << next_kmer[-1] # add only last base
    end

    hash[kmer].delete next_kmer

    hash.delete(kmer) if hash[kmer].empty?

  else
    # check if all edges have been hit
    arr = hash.select { |kmer, kmers| kmers.count > 0 }

    if arr.empty?
      contigs[num] = tmp_walk
      num += 1
      break
    else
      contigs[num] = tmp_walk
      num += 1

      if branch_point_kmers.any? { |kmer, info| !info[:end_of].empty? }
        next_kmer = branch_point_kmers.first.first
        branch_point_kmers[next_kmer][:start_of] << num
      else

        next_kmer = in_out.select { |kmer, connections| connections[:in].zero? }
        if next_kmer.nil? || next_kmer.empty?
          next_kmer = arr.first.first
        else
          next_kmer = next_kmer.first.first
        end
      end

      tmp_walk = String.new(next_kmer)
    end
  end

  kmer = next_kmer
end

inbase = File.basename(ARGV.first, File.extname(ARGV.first))
outbase = File.join(File.dirname(ARGV.first), "shaq.#{inbase}")
contigs_f = "#{outbase}.contigs.fa"
graph_f = "#{outbase}.contig_graph.txt"

File.open(contigs_f, "w") do |f|
  contigs.each do |num, contig|
    f.puts ">contig_#{num}\n#{contig}"
  end
end

File.open(graph_f, "w") do |f|
  f.puts %w[C1 C2 BranchingKmer].join "\t"
  branch_point_kmers.each do |kmer, info|
    info[:end_of].each do |end_of|
      info[:start_of].each do |start_of|
        f.puts ["contig_#{end_of}", "contig_#{start_of}", kmer].join "\t"
      end
    end
  end
end
