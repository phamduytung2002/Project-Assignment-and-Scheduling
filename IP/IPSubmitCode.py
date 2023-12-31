from ortools.linear_solver import pywraplp

#read input from file
#Integer programming model
def ILP():
    n, q = map(int, input().split())
    worker_list =[[] for _ in range(n)]
    d = []
    pres = []
    starts = []
    for i in range(q):
        pres.append(tuple(map(int, input().split())))
    d = list(map(int, input().split()))
    m = int(input())
    starts = list(map(int, input().split()))
    K = int(input())
    key = 1
    dat = []
    cost_list = [[999 for _ in range(m)] for _ in range(n)]
    for i in range(K):
        line = input()
        parts = line.split()
        if len(parts) == 3:
            i, j, value = map(int, parts)
            worker_list[i - 1].append(j-1)
            cost_list[(i-1)][(j-1)] = value
    dct = {}
    for i in pres:
        if i[0]-1 not in dct:
            dct[i[0]-1] = [i[1]-1]
        else:
            dct[i[0]-1].append(i[1]-1)

    solver = pywraplp.Solver.CreateSolver('SCIP')
    INF = solver.infinity()
    large = max(sum(d) + max(d)+max(starts), n + 1)
    # X = m*n
    X = [[None]*n for i in range(m)]
    
    # X[i][j] team i vs task j
    for i in range(m):
        for j in range(n):
            X[i][j] = solver.IntVar(0, 1, 'X[{}, {}]'.format(i, j))

    for i in range(n):
        # Sua lai dieu kien
        solver.Add(sum([X[j][i] for j in range(m)]) <= 1)

    for i in range(n):
        not_in = set(list(range(m))) - set(worker_list[i])
        for j in not_in:
            solver.Add(X[j][i] == 0)

    assign = [solver.IntVar(0, m-1, 'A[{}]'.format(i)) for i in range(n)]

    for i in range(m):
        for j in range(n):
            solver.Add(large * (-X[i][j]+1) + i >= assign[j])
            solver.Add(-large * (-X[i][j]+1) + i <= assign[j])

    times = [solver.IntVar(0, large, 'times[{}]'.format(i)) for i in range(n)]

    for first in dct:
        for second in dct[first]:
            solver.Add(times[first] + d[first] <= times[second])

    for i in range(m):
        for j in range(n):
            for k in range(j+1, n):
                if j != k:
                    t = solver.IntVar(0, 1, 'T[{}, {}, {}]'.format(i, j, k))
                    solver.Add((1 - X[i][j])*large + (1 - X[i][k])*large + times[j] + t*large >= times[k] + d[k])
                    solver.Add((1 - X[i][j])*large + (1 - X[i][k])*large + times[k] + (1-t)*large >= times[j] + d[j])
            solver.Add((1-X[i][j])*large + times[j] >= starts[i])

    costs = solver.IntVar(0, INF, 'cost')
    sum_ = 0
    for i in range(m):
        for j in range(n):
            sum_ += X[i][j] * cost_list[j][i]
    solver.Add(costs == sum_)
    number_task = solver.IntVar(1, n, 'n_t')
    task = 0
    for i in range(m):
        for  j in range(n):
            task += X[i][j]
    solver.Add(number_task == task)
    Z = solver.IntVar(0, INF, 'z')
    for i in range(n):
        solver.Add(Z >= times[i] + d[i])
    solver.Maximize(number_task)
    status = solver.Solve()
    if status == pywraplp.Solver.OPTIMAL:
        solver.Add(number_task == number_task.solution_value())
        solver.Minimize(Z)
        status = solver.Solve()
        if status == pywraplp.Solver.OPTIMAL:
            solver.Add(Z == Z.solution_value())
            solver.Minimize(costs) 
            status = solver.Solve()
            if status == 0:
                print(int(number_task.solution_value()))
                # for i in range(n):
                #     if X[int(assign[i].solution_value())][i].solution_value() ==1:
                #         print('Time for part {}:'.format(i+1), times[i].solution_value())
                for i in range(n):
                    if X[int(assign[i].solution_value())][i].solution_value() ==1:
                        print(i+1, int(assign[i].solution_value()+1), int(times[i].solution_value()))
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
if __name__ == '__main__':
    ILP()