"""
   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
"""

#
# Problem from the film "Good Will Hunting"
#
#   "Draw all homomorphically irreducible trees for n = 10"
#

from docplex.cp.model import CpoModel, CpoSolver
import sys
import subprocess
import os


def gwh(n, display_screen=True, display_graphs=True):
    # Variables
    mdl = CpoModel()
    x = [[0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            x[i][j] = x[j][i] = mdl.integer_var(0, 1, 'X({}, {})'.format(i, j))

    # Nodes must be connected (degree != 0) and not contractable (degree != 2)
    deg = [mdl.sum(x[i]) for i in range(n)]
    if n > 1:
        for i in range(n):
            mdl.add(mdl.forbidden_assignments(deg[i], [0, 2]))

    # Avoid cycles. My rank is r.  One neighbor (my parent) should have
    # rank r-1.  Other neighbors (children) should have rank r+1.
    rank = [mdl.integer_var(0, n - 1, 'R({})'.format(i)) for i in range(n)]
    mdl.add(rank[0] == 0)
    for i in range(1, n):
        num_parents = 0
        for j in range(n):
            mdl.add(mdl.if_then(x[i][j] == 1, mdl.abs(rank[i] - rank[j]) == 1))
            num_parents += (x[i][j] == 1) & (rank[j] == rank[i] - 1)
        mdl.add(num_parents == 1)

    # Symmetry breaking.  Order pairs of nodes according to degree and then
    # connection signature (not counting the pair of nodes under consideration)
    # This is enough to produce the correct number of 10 solutions for gwh(10),
    # 14 for gwh(11). For gwh(12), 28 solutions are found (2 duplicate trees
    # as there are actually 26).
    for i in range(1, n):
        a1 = [deg[i-1]] + x[i-1][:i-1]   + x[i-1][i+1:]
        a2 = [deg[i  ]] + x[i  ][:i-1]   + x[i  ][i+1:]
        mdl.add(mdl.lexicographic(a2, a1))

    # Open .dot file
    if display_graphs:
        dotfile = 'gwh_{}.dot'.format(n)
        dot = open(dotfile, "w")
        print('graph {', file=dot)

    # Iterate over all the solutions, and display
    s = 0
    for res in CpoSolver(mdl, LogVerbosity='Quiet'):
        s += 1
        if display_screen:
            sig = [sum(res[x[i][j]] for j in range(n) if i != j)
                   for i in range(n)]
            print('Graph {}, degrees = {}'.format(s, sig))
        for i in range(n):
            line = '  ' * i + '.'
            for j in range(i + 1, n):
                line += ' {}'.format(res[x[i][j]])
                if display_graphs and res[x[i][j]]:
                    ni, nj = chr(ord('A') + i), chr(ord('A') + j)
                    print('"{0}{1}"--"{0}{2}";'.format(s, ni, nj), file=dot)
            if display_screen:
                print(line)

    # Finish up
    print('{} solutions found'.format(s))
    if display_graphs:
        print('}', file=dot)
        dot.close()
        with open('gwh_{}.pdf'.format(n), "w") as pdf:
            subprocess.run(['neato', '-Tpdf', 'gwh_{}.dot'.format(n)], stdout=pdf)
        os.unlink(dotfile)


if __name__ == '__main__':
    n = 10  # Problem size from the film
    pdf_out = True

    if len(sys.argv) > 1:
        n = int(sys.argv[1])
    try:
        subprocess.run(['neato', '-h'])
    except FileNotFoundError as e:
        pdf_out = False

    gwh(n, display_graphs=pdf_out)
    if pdf_out:
        print("GraphViz rendered trees in gwh_{}.pdf".format(n))
    else:
        print("GraphViz rendering commands not found.  No PDF generated")
