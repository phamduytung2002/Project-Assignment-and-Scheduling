#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>

using namespace std;

const int MAXN = 1000;
const int MAXM = 500;
const int MAXK = 1000000;

int N, Q, M, K;
int d[MAXN+1];
int s[MAXM+1];
int c[MAXN+1][MAXM+1];
vector<int> adj[MAXN+1];
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
            return scheduled.size() > other.scheduled.size();
        }
        if (earliest_completion_time != other.earliest_completion_time) {
            return earliest_completion_time > other.earliest_completion_time;
        }
        return cost > other.cost;
    }
};

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
    while (!pq.empty()) {
        State state = pq.top();
        pq.pop();
        if (state.scheduled.size() == N) {
            cout << state.scheduled.size() << endl;
            for (int i = 0; i < N; i++) {
                cout << state.scheduled[i] << " " << state.assigned[i] << " " << state.start_times[i] << endl;
            }
            return;
        }
        for (int i = 1; i <= N; i++) {
            if (find(state.scheduled.begin(), state.scheduled.end(), i) == state.scheduled.end()) {
                for (int j = 1; j <= M; j++) {
                    if (c[i][j] > 0) {
                        State child = state;
                        child.task = i;
                        child.team = j;
                        // Compute the start time for the task considering its dependencies and the availability of the team
                        int start_time = max(child.time, child.team_availability[j]);
                        for (int pred : adj[i]) {
                            auto it = find(child.scheduled.begin(), child.scheduled.end(), pred);
                            if (it != child.scheduled.end()) {
                                int pred_index = it - child.scheduled.begin();
                                // Use the finish time of the predecessor task
                                start_time = max(start_time, child.start_times[pred_index] + d[child.scheduled[pred_index]]);
                            }
                        }
                        child.start_times.push_back(start_time);
                        // Update the availability time of the team
                        child.team_availability[j] = start_time + d[i];
                        // Update the current time in the state
                        child.time = max(child.time, start_time + d[i]);
                        child.scheduled.push_back(i);
                        child.assigned.push_back(j);
                        pq.push(child);
                    }
                }
            }
        }
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
    branch_and_bound();
    return 0;
}
