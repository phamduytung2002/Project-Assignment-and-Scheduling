from ortools.sat.python import cp_model

def CP():
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

    model = cp_model.CpModel()

    summ = 0
    for i in cost_list:
        summ += sum(i)
    summ += sum(starts) 
    summ += sum(d)

    INF = summ
    X = [[None]*n for i in range(m)]

    for i in range(m):
        for j in range(n):
            X[i][j] = model.NewIntVar(0, 1, 'X[{}, {}]'.format(i, j))

    for i in range(n):
        model.Add(sum([X[j][i] for j in range(m)]) <= 1) #constrain each part of problem only has 1 worker

    for i in range(n):
        not_in = set(list(range(m))) - set(worker_list[i])
        for j in not_in:
            model.Add(X[j][i] == 0) #constrain the worker i can not do the part of problem j

    assign = [model.NewIntVar(0, m-1, 'A[{}]'.format(i)) for i in range(n)]
    for i in range(m):
        for j in range(n):
            model.Add(assign[j] == i).OnlyEnforceIf(X[i][j])

    times = [model.NewIntVar(0, sum(d) + max(starts), 'times[{}]'.format(i)) for i in range(n)]

    for first in dct:
        for second in dct[first]:
            model.Add(times[first] + d[first] <= times[second]) #constrain order of part of problem


    for i in range(m):
        for j in range(n):
            model.Add(times[j] >= starts[i]).OnlyEnforceIf(X[i][j])
            for k in range(j+1, n):
                b = model.NewBoolVar('b')
                model.Add(times[j] + d[j] <= times[k] ).OnlyEnforceIf(X[i][j]).OnlyEnforceIf(X[i][k]).OnlyEnforceIf(b)
                model.Add( times[k] + d[k] <= times[j] ).OnlyEnforceIf(X[i][j]).OnlyEnforceIf(X[i][k]).OnlyEnforceIf(b.Not())
    costs = model.NewIntVar(0, INF, 'cost')
    sum_ = 0
    for i in range(m):
        for j in range(n):
            sum_ += X[i][j] * cost_list[j][i]
    model.Add(costs == sum_)
    
    number_task = model.NewIntVar(0, n, 'n_t')
    task = 0
    for i in range(m):
        for  j in range(n):
            task += X[i][j]
    model.Add(number_task == task)
            
    Z = model.NewIntVar(0, INF, 'z')
    for i in range(n):
        model.Add(Z >= times[i] + d[i])
        
    model.Maximize(number_task)
    solver = cp_model.CpSolver()
    status = solver.Solve(model)
    number_task_max = solver.Value(number_task)
    if status == cp_model.OPTIMAL:
        model.Add(number_task == number_task_max)
        model.Minimize(Z)
        status  = solver.Solve(model)
        Z_min= solver.Value(Z)
        if status == cp_model.OPTIMAL:
            model.Add(Z == Z_min)
            model.Minimize(costs)
            status  = solver.Solve(model)
            if status == 4:
                print(solver.Value(number_task))
                # print('      Minimum time:', solver.Value(Z))
                # print('      Minimum cost:', solver.Value(costs))
                # for i in range(n):
                #     print('Starting time for part {}:'.format(i+1), solver.Value(times[i]))

                for i in range(n):
                    # if X[int(solver.Value(assign[i]))][i] == 1:
                    print(i+1, solver.Value(assign[i])+1, solver.Value(times[i]))
            else:
                if status == 3:
                    print('INFEASIBLE')
                elif status == 2:
                    print('FEASIBLE')
                elif status == 1:
                    print('MODEL_INVALID')
                elif status == 0:
                    print('UNKNOWN')
                else:
                    print('ERROR')
        else:
                if status == 3:
                    print('INFEASIBLE')
                elif status == 2:
                    print('FEASIBLE')
                elif status == 1:
                    print('MODEL_INVALID')
                elif status == 0:
                    print('UNKNOWN')
                else:
                    print('ERROR')
if __name__ == '__main__':
    CP()
