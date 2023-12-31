import numpy as np
import math
import time

# n: Số lượng công nhân
# q: Số lượng công việc
# worker_list[i]: Danh sách công việc cho mỗi công nhân
# pres: Phụ thuộc giữa các công việc với các dòng tiếp theo
# d: Thời gian thực hiện của từng công việc
# m: Số lượng nhóm công nhân
# starts[j]: đọc thời gian bắt đầu công việc của từng đội
# K: Đọc số lượng công việc chính
# cost_list: Khởi tạo ma trận chi phí với giá trị ban đầu là 999 
def input(filename):
    with open(filename, 'r') as f:
        n, q = map(int, f.readline().split())               
        worker_list = [[] for _ in range(n)]
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
            lines = f.readlines()
            for line in lines:
                parts = line.split()
                if len(parts) == 3:
                    i, j, value = map(int, parts)
                    worker_list[i-1].append(j-1)  # Cập nhật danh sách công việc cho từng đội
                    cost_list[(i-1)][(j-1)] = value  # Cập nhật chi phí của từng công việc
    return n, m, worker_list, d, K, cost_list, starts, pres

def greedy_algorithm(filename):
    starting = time.time()
    # Đọc dữ liệu từ file
    n, m, worker_list, d, K, cost_list, starts, pres = input(filename)
    
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
    print(task_list)
    for i in range(n):
        if i+1 not in task_list:
            task_list.append(i+1)
    print(task_list)
    loop = 0
    exit = 0
    print("Assigning tasks to workers...")
    minimum_cost = 0;
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
    if exit != 1:
        print('------------------------------------------')
        print('OPTIMAL SOLUTION FOUND!  ')
        for i in range(n):
            print('Task ' + str(i+1) + ' is performed by group ' + str(solution[i][0] + 1) + ' from time ' + str(solution[i][1]) + ' to ' + str(solution[i][2]) + ' with cost = ' + str(solution[i][3]))
        print('Min time to finish: ' + str(final_time))
        print('Min cost: ' + str(minimum_cost))


    ending = time.time()
    total_time = ending - starting
    print('------------------------------------------')
    print('Execution time: ' + str(total_time))

