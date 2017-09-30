#!/usr/bin/ruby -w

##############################################
#
# preprocessing the migra1 .dat files
#
##############################################

y = "-1"

ARGV.each do |i|
  puts i
end 


if ARGV[0] == "-y"
  y = ARGV[1]
  f_name = ARGV[2]
else 
  f_name = ARGV[0]
end

max_exp = 2 # maximal experience
year = -1
n_year = 0

sched = Hash.new
gyak = Hash.new
surv = Hash.new
dat = File.open(f_name,"r")
dat.each_line do |line|
  v = line.split
  if v[0] == y || y == "-1"
    if year != v[0] 
      year = v[0]
      n_year += 1
    end 
    if sched[v[2]] == nil 
      sched[v[2]] = ""
      gyak[v[2]] = 0
      surv[v[2]] = 0
    end
    if v[1] == 0 && v[8] == max_exp
      gyak[v[2]] += 1
    end
    if v[1] == "0" 
      sched[v[2]] = sched[v[2]] + "|W,#{v[1]},#{v[10]},#{v[6]},#{v[7]},#{v[11]}"
    elsif v[1] == "26" 
      sched[v[2]] = sched[v[2]] + "|U,#{v[1]},#{v[10]},#{v[6]},#{v[7]},#{v[11]}"
    end
    if v[4] =~ /^S$/ 
      sched[v[2]] = sched[v[2]] + "|M,#{v[1]},#{v[10]},#{v[6]},#{v[7]},#{v[11]}"
    elsif v[4] =~ /^E$/ 
      sched[v[2]] = sched[v[2]] + "|e,#{v[1]},#{v[10]},#{v[6]},#{v[7]},#{v[11]}"
    end
    if v[5] =~ /^No./ 
      sched[v[2]] = sched[v[2]] + "|N,#{v[1]},#{v[10]},#{v[6]},#{v[7]},#{v[11]},#{v[12]},#{v[13]}"
    elsif v[5] =~ /^So./ 
      sched[v[2]] = sched[v[2]] + "|S,#{v[1]},#{v[10]},#{v[6]},#{v[7]},#{v[11]},#{v[12]},#{v[13]}"
    end
    if v[5] =~ /^S$/ 
      sched[v[2]] = sched[v[2]] + "|B,#{v[1]},#{v[10]},#{v[6]},#{v[7]},#{v[11]}"
    elsif v[5] =~ /^A$/ 
      sched[v[2]] = sched[v[2]] + "|a,#{v[1]},#{v[10]},#{v[6]},#{v[7]},#{v[11]}"
    elsif v[5] =~ /^AT$/ 
      sched[v[2]] = sched[v[2]] + "|d,#{v[1]},#{v[10]},#{v[6]},#{v[7]},#{v[11]}"
    elsif v[5] =~ /^B$/
      sched[v[2]] = sched[v[2]] + "|T,#{v[1]},#{v[10]},#{v[6]},#{v[7]},#{v[11]}"
    elsif v[5] =~ /.B$/
      sched[v[2]] = sched[v[2]] + "|T,#{v[1]},#{v[14]},#{v[6]},#{v[7]},#{v[11]}"
    end
    if v[1] == 51 && gyak[v[2]] == n_year
      surv[v[2]] += 1
    end 
  end
end
dat.close

for bird in surv do
  if surv[bird] == n_year
    if sched[bird] == "" 
      puts "#{bird} n"
    else 
      puts "#{bird} #{sched[bird]}"
    end
  end
end
