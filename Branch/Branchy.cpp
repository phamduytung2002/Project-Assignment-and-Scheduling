#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>
#include <climits>
#include <tuple>

using namespace std;

const int MAXN = 1000;
const int MAXM = 500;
const int MAXK = 1000000;

int N, Q, M, K;
vector<int> d(MAXN+1), s(MAXM+1);
vector<vector<int>> adj(MAXN+1), c(MAXN+1, vector<int>(MAXM+1));
vector<int> team_availability(MAXM+1, 0);
struct State {
    int task;
    int team;
    int time;
    vector<int> scheduled;
    vector<int> assigned;
    vector<int> start_times;
    vector<int> team_availability;
    int cost;
    int earliest_completion_time;

    bool operator<(const State& other) const {
        if (scheduled.size() != other.scheduled.size()) {
            return scheduled.size() < other.scheduled.size();
        }
        if (earliest_completion_time != other.earliest_completion_time) {
            return earliest_completion_time < other.earliest_completion_time;
        }
        return cost < other.cost;
    }
};

vector<int> visited(MAXN+1, 0);
vector<int> recStack(MAXN+1, 0);

bool isCyclicUtil(int i, vector<int>& tasksToRemove) {
    if (!visited[i]) {
        visited[i] = 1;
        recStack[i] = 1;

        for (int j : adj[i]) {
            if (!visited[j] && isCyclicUtil(j, tasksToRemove))
                return true;
            else if (recStack[j]) {
                tasksToRemove.push_back(j);
                return true;
            }
        }
    }
    recStack[i] = 0;
    return false;
}

void removeCyclicTasks() {
    vector<int> tasksToRemove;
    for (int i = 1; i <= N; i++)
        if (isCyclicUtil(i, tasksToRemove))
            for (int task : tasksToRemove)
                adj[task].clear();  // Remove all dependencies of the task

    // Remove tasks from the list of tasks
    for (int task : tasksToRemove)
        d[task] = -1;  // Mark the task as removed
}

int compute_earliest_completion_time(const State& state) {
    vector<int> earliest_start_time(N+1, 0);
    for (int task : state.scheduled) {
        int earliest_start = 0;
        for (int pred : adj[task]) {
            earliest_start = max(earliest_start, earliest_start_time[pred] + d[pred]);
        }
        earliest_start_time[task] = earliest_start;
    }
    int critical_path = 0;
    for (int task : state.scheduled) {
        critical_path = max(critical_path, earliest_start_time[task] + d[task]);
    }
    return critical_path + state.time;
}

void branch_and_bound() {
    priority_queue<State> pq;
    pq.push({0, 0, 0, {}, {}, {}, vector<int>(MAXM+1, 0), 0, compute_earliest_completion_time({0, 0, 0, {}, {}, {}, vector<int>(MAXM+1, 0), 0, 0})});
    int min_cost = INT_MAX;
    State best_state;
    while (!pq.empty()) {
        State state = pq.top();
        pq.pop();
        if (state.cost >= min_cost) continue;
        if (state.scheduled.size() == N) {
            min_cost = state.cost;
            best_state = state;
            continue;
        }
        for (int i = 1; i <= N; i++) {
            if (find(state.scheduled.begin(), state.scheduled.end(), i) == state.scheduled.end()) {
                bool all_predecessors_scheduled = true;
                for (int pred : adj[i]) {
                    if (find(state.scheduled.begin(), state.scheduled.end(), pred) == state.scheduled.end()) {
                        all_predecessors_scheduled = false;
                        break;
                    }
                }
                if (!all_predecessors_scheduled) continue;
                for (int j = 1; j <= M; j++) {
                    if (c[i][j] > 0) {
                        State child = state;
                        child.task = i;
                        child.team = j;
                        int start_time = max({child.team_availability[j], s[j]});
                        for (int pred : adj[i]) {
                            auto it = find(child.scheduled.begin(), child.scheduled.end(), pred);
                            if (it != child.scheduled.end()) {
                                int pred_index = it - child.scheduled.begin();
                                start_time = max(start_time, child.start_times[pred_index] + d[child.scheduled[pred_index]]);
                            }
                        }
                        child.start_times.push_back(start_time);
                        child.team_availability[j] = start_time + d[i];
                        child.time = min(child.time, start_time + d[i]);
                        child.scheduled.push_back(i);
                        child.assigned.push_back(j);
                        child.cost += c[i][j];
                        child.earliest_completion_time = compute_earliest_completion_time(child);
                        if (child.cost + child.earliest_completion_time >= min_cost) continue;
                        pq.push(child);
                    }
                }
            }
        }
    }
    vector<tuple<int, int, int>> output;
    for (int i = 0; i < N; i++) {
        output.push_back(make_tuple(best_state.scheduled[i], best_state.assigned[i], best_state.start_times[i]));
    }
    sort(output.begin(), output.end());
    cout << output.size() << endl;
    for (auto& out : output) {
        cout << get<0>(out) << " " << get<1>(out) << " " << get<2>(out) << endl;
    }
}

int main() {
    cin >> N >> Q;
    for (int i = 1; i <= Q; i++) {
        int a, b;
        cin >> a >> b;
        adj[b].push_back(a);
    }
    for (int i = 1; i <= N; i++) {
        cin >> d[i];
    }
    cin >> M;
    for (int i = 1; i <= M; i++) {
        cin >> s[i];
    }
    cin >> K;
    for (int i = 1; i <= K; i++) {
        int a, b, cost;
        cin >> a >> b >> cost;
        c[a][b] = cost;
    }
    removeCyclicTasks();
    branch_and_bound();
    return 0;
}
