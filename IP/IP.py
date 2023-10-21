from ortools.linear_solver import pywraplp

#read input from file
def input(filename):
    with open(filename, 'r') as f:
        n, q = map(int, f.readline().split())
        worker_list =[]
        d = []
        pres = []
        for i in range(q):
            pres.append(tuple(map(int, f.readline().split())))
        d = list(map(int, f.readline().split()))
        m = int(f.readline())
        starts = list(map(int, f.readline().split()))
        K = int(f.readline())
        key = 1
        dat = []
        cost_list = [[999 for _ in range(m)] for _ in range(n)]
        for i in range(K):
            read_line = list(map(int, f.readline().split()))
            if read_line [0] == key and i < (K-1):
                dat.append(read_line[1])
            elif (key != read_line[0]) or (i == (K-1)):
                if key == read_line[0]:
                    dat.append(read_line[1]-1)
                worker_list.append(dat)
                dat =[]
                key += 1
                dat.append(read_line[1]-1)
            cost_list[read_line[0]-1][read_line[1]-1] = read_line[2]
            
    return n, m, worker_list, d, K, cost_list, starts, pres

#Integer programming model
def ILP(filename):
    n, m, worker_list, d, K, c, starts, pres = input(filename)
    dct = {}
    for i in pres:
        if i[0]-1 not in dct:
            dct[i[0]-1] = [i[1]-1]
        else:
            dct[i[0]-1].append(i[1]-1)

    solver = pywraplp.Solver.CreateSolver('CBC')
    INF = solver.infinity()
    large = max(sum(d) + max(d)+max(starts), n + 1)

    X = [[None]*n for i in range(m)]

    for i in range(m):
        for j in range(n):
            X[i][j] = solver.IntVar(0, 1, 'X[{}, {}]'.format(i, j))

    for i in range(n):
        solver.Add(sum([X[j][i] for j in range(m)]) == 1)

    for i in range(n):
        not_in = set(list(range(m))) - set(worker_list[i])
        for j in not_in:
            solver.Add(X[j][i] == 0)

    assign = [solver.IntVar(0, m-1, 'A[{}]'.format(i)) for i in range(n)]

    for i in range(m):
        for j in range(n):
            solver.Add(large * (-X[i][j]+1) + i >= assign[j])
            solver.Add(-large * (-X[i][j]+1) + i <= assign[j])

    times = [solver.IntVar(0, sum(d), 'times[{}]'.format(i)) for i in range(n)]

    for first in dct:
        for second in dct[first]:
            solver.Add(times[first] + d[first] <= times[second])

    for i in range(m):
        for j in range(n):
            for k in range(j+1, n):
                if j != k:
                    t = solver.IntVar(0, 1, 'T[{}, {}, {}]'.format(i, j, k))
                    solver.Add((1 - X[i][j])*large + (1 - X[i][k])*large + times[j] + (1-t)*large >= times[k] + d[k])
                    solver.Add((1 - X[i][j])*large + (1 - X[i][k])*large + times[k] + t*large >= times[j] + d[j])
            solver.Add((1-X[i][j])*large + times[j] >= starts[i])


    costs = [solver.IntVar(0, INF, 'cost[{}]'.format(i)) for i in range(m)]
    for i in range(m):
        sum_ = 0
        for j in range(n):
            sum_ += X[i][j] * c[j][i]
        solver.Add(costs[i] == sum_)
    number_task = solver.IntVar(0, n, 'n_t')
    solver.Add(number_task <= n)
    Z = solver.IntVar(0, INF, 'z')
    max_cost = solver.IntVar(0, INF, 'maxcost')
    for i in range(n):
        solver.Add(Z >= times[i] + d[i])
    for i in range(m):
        solver.Add(max_cost >= costs[i])
    solver.Maximize(number_task)
    status = solver.Solve()
    if status == pywraplp.Solver.OPTIMAL:
        solver.Add(number_task == number_task.solution_value())
        solver.Minimize(Z)
        status = solver.Solve()
        if status == pywraplp.Solver.OPTIMAL:
            solver.Add(Z == Z.solution_value())
            solver.Minimize(max_cost)
            status = solver.Solve()
            if status == 0:
                
                print('Optimal result: ')
                print('   Maximum task:', number_task.solution_value())
                print('   Maximum cost:', max_cost.solution_value())
                print('   Minimum time:', Z.solution_value())
                for i in range(n):
                    print('Time for part {}:'.format(i+1), times[i].solution_value())

                for i in range(n):
                    print('The worker group of part {}:'.format(i+1), assign[i].solution_value()+1)
            else:
                if status == 2:
                    print('INFEASIBLE')
                elif status == 3:
                    print('Unbound')
                else:
                    print('ERROR')
        else:
            if status == 2:
                print('INFEASIBLE')
            elif status == 3:
                print('Unbound')
            else:
                print('ERROR')   