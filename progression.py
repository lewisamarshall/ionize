import ionize

sol = ionize.Solution(['tris', 'acetic acid'], [0.01, 0.005])

new_sol = sol.displace('tris', 'bis-tris')
# print(sol.concentrations, sol.ions)
reversion = new_sol.displace('bis-tris', 'tris')
# print(new_sol.displace('bis-tris', 'tris').ions)

for s in (sol, new_sol, reversion):
    print('alberty: {}'.format(s.alberty()))
    print('jovin: {}'.format(s.jovin()))
    print('concentrations: {}'.format(s.concentrations))
