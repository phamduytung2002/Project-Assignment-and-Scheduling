import numpy as np
import math
import time
def greedy_algorithm():
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

    starting = time.time()
    task_list = []
    cost = 0
    final_time = 0
    timefinish = {}
    solution = {}
    cost = [0] * m

    # Xác định danh sách công việc cần thực hiện trước
    for (u, v) in pres:
        if u not in task_list:
            task_list.append(u)
    for i in range(n):
        if i+1 not in task_list:
            task_list.append(i+1)
    loop = 0
    exit = 0
    minimum_cost = 0;
    task_count = 0
    while len(task_list):
        consider_task = task_list[0]
        count = 0
        if loop == n+1:
            print("No solution")
            exit = 1
            break

        # Nếu công việc có công việc tiên quyết, đưa nó xuống cuối danh sách
        for u in range(n):
            if (u, consider_task) in pres and u in task_list:
                task_list.pop(0)
                task_list.append(consider_task)
                count += 1
                loop += 1
                break;
        if count == 1:
            continue
        loop = 0
        consider_task = task_list[0]
        min_start = starts[0]
        min_worker = 0
        
        # Tối ưu thời gian
        for worker in worker_list[consider_task-1]:
            if starts[worker] < min_start:
                min_start = starts[worker]
                min_worker = worker
        
        # Tính chi phí nếu gán công việc này cho các nhóm
        cost_list_extend = [cost[i] + cost_list[consider_task-1][i] for i in range(m)]
        min_cur_cost = cost_list_extend[min_worker]
        # Tối ưu chi phí
        for worker in worker_list[consider_task-1]:
            if starts[worker] == min_start and cost_list_extend[worker] < min_cur_cost:
                min_worker = worker
                min_cur_cost = cost_list_extend[worker]

        # So sánh thời gian công việc này với thời gian hoàn thành của công việc tiên quyết
        for (u, v) in pres:
            if consider_task == v and u in timefinish.keys():
                if min_start < timefinish[u]:
                    min_start = timefinish[u]

        cost[min_worker] += cost_list[consider_task-1][min_worker]
        starts[min_worker] = min_start + d[consider_task-1]
        minimum_cost += cost_list[consider_task-1][min_worker]
        solution[consider_task-1] = [min_worker, min_start, starts[min_worker], cost_list[consider_task-1][min_worker]]

        if starts[min_worker] > final_time:
            final_time = starts[min_worker]
        timefinish[consider_task] = starts[min_worker]
        task_list.pop(0)
        task_count += 1
    if exit != 1:
        print(task_count)
        for i in range(n):
            print(f"{i+1} {solution[i][0]+1} {solution[i][1]}")

if __name__ == '__main__':
    greedy_algorithm()