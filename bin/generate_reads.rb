readsize = 25
numreads = 100
dna = "i_love_to_eat_apple_pie_and_love_to_eat_other_fun_things"

start_range = (0..dna.length-readsize).to_a
reads = numreads.times.each_with_index do |n|
  start = start_range.sample

  puts [">read_#{n}", dna[start, readsize]].join "\n"
end
