import sys

testlog_name = sys.argv[1].split(".")[0]

ls = [l.strip().split(",") for l in open(f"{testlog_name}.csv").readlines()]

dic = { (int(l[0]), int(l[1])): list(map(float, l[2:])) for l in ls }
ds = list(sorted(set([int(l[0]) for l in ls])))

ks = [3, 5, 10, 15, 20]

def write_log(ls, header, name):
  f = open(f'{testlog_name}_{name}.csv', "w")
  assert len(header) == len(ls[0])
  print(*header, sep=",", file=f)
  for l in ls:
    print(*l, sep=",", file=f)

# proposed_precalc
proposed_precalc = [
  [d] + 
  [dic[d, k][1] for k in ks]
  for d in ds
]
write_log(proposed_precalc, ["d", *[f'precalc_{k}' for k in ks]], "precalc")

# half_gcd elapsed / proposed_total
total_compare = [
  [d] + 
  [sum([dic[d, k][0] for k in ks]) / 5] +
  [dic[d, k][4] for k in ks]
  for d in ds
]
write_log(total_compare, ["d", "half_gcd", *[f'proposed_{k}' for k in ks]], "total_compare")

# each details
for k in ks:
  each_compair = [[d, dic[d, k][2], dic[d, k][3], dic[d, k][4]] for d in ds]
  write_log(each_compair, ["d", "compute_Mv", "compute_gcd", "total"], f"details_{k}")
