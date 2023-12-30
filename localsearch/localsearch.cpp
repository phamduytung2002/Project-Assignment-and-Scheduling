// #define NDEBUG
#include <bits/stdc++.h>
using namespace std;

const int VERY_BIG_NUMBER = 1e7;
const int PENALTY_PER_FAULT = 1e3;
const double MAXTIME = 100; // second
const int MAXITER = 1000;
const int maxn = 1005;
const int maxm = 505;
const int maxq = 1755;
const int MAX_PATIENT = 20;

bool must_before[maxn * maxn];

class Problem {
  public:
    Problem() {}

    // given
    int N;                  // number of tasks
    int M;                  // number of teams
    int s[maxm];            // team j start at s[j]
    int K;                  // number of (team, task) pairs
    map<int, int> c[maxm];  // team j do i cost c[i][j]
    int d[maxm];            // task i has duration d[i]
    int Q;                  // number of precedence constraints
    pair<int, int> p[maxq]; // task [i].first must be done before p[i].second

    // addition
    vector<int>
        can_do_by_team[maxm]; // team j can do tasks in can_do_by_team[j]
    vector<int> succ[maxn]; // task i must be done before every tasks in succ[i]
    // bool must_before[maxnm * maxnm]; // must_before[i][j] = true if task i
    // must
    //                                  // be done before task j, or
    //                                  // i is an ancestor of j
    vector<int>
        tasks_toposorted; // an order of tasks can be done by a single team

    bool get_must_before(int row, int col) {
        return must_before[row * maxn + col];
    }

    void set_must_before(int row, int col, bool val) {
        must_before[row * maxn + col] = val;
    }

    void input() {
        cin >> N >> Q;
        for (int i = 1; i <= Q; ++i) {
            // input precedence constraints
            int u, v;
            cin >> u >> v;
            p[i] = make_pair(u, v);
            succ[u].push_back(v);
        }
        for (int i = 1; i <= N; ++i)
            cin >> d[i];
        cin >> M;
        for (int i = 1; i <= M; ++i)
            cin >> s[i];
        cin >> K;

        for (int i = 1; i <= K; ++i) {
            int u, v, w;
            cin >> u >> v >> w;
            c[u][v] = w;
            can_do_by_team[v].push_back(u);
        }
    }
} prob;

class Solution {
  private:
    void check_all_descendant_cantdo(int u) {
        for (int i = 0; i < prob.succ[u].size(); ++i) {
            int v = prob.succ[u][i];
            if (cando[v]) {
                cando[v] = false;
                check_all_descendant_cantdo(v);
            }
        }
    }

    void dfs(int u, bool visited[], bool rec[]) {
        // cando[u] = false;
        for (int i = 0; i < prob.succ[u].size(); ++i) {
            int v = prob.succ[u][i];
            if (cando[v]) {
                if (!visited[v]) {
                    visited[v] = true;
                    rec[v] = true;
                    dfs(v, visited, rec);
                    rec[v] = false;
                } else if (rec[v]) {
                    cando[v] = false;
                    check_all_descendant_cantdo(v);
                }
            }
        }
    }

    void update_cando() {
        // cando[i] = false if task i is included in a cycle
        // after this we can ignore every task i that cando[i] = false
        for (int i = 1; i <= prob.N; ++i)
            cando[i] = true;
        bool visited[prob.N + 1];
        memset(visited, 0, sizeof(visited));
        bool rec[prob.N + 1];
        memset(rec, 0, sizeof(rec));
        for (int i = 1; i <= prob.N; ++i) {
            if (!visited[i])
                dfs(1, visited, rec);
        }
        for (int i = 1; i <= prob.N; ++i) {
            R += cando[i];
        }
    }

    void dfs_update_must_before(int u) {
        for (int i = 0; i < prob.succ[u].size(); ++i) {
            int v = prob.succ[u][i];
            if (cando[v]) {
                dfs_update_must_before(v);
            }
            for (int j = 1; j <= prob.N; ++j) {
                if (prob.get_must_before(v, j)) {
                    prob.set_must_before(u, j, true);
                }
            }
            prob.set_must_before(u, v, true);
        }
    }

    void update_must_before() {
        for (int i = 1; i <= prob.N; ++i) {
            if (cando[i])
                dfs_update_must_before(i);
        }
    }

    void dfs_toposort(int u, bool visited[]) {
        visited[u] = true;
        for (int i = 0; i < prob.succ[u].size(); ++i) {
            int v = prob.succ[u][i];
            if (!visited[v] && cando[v]) {
                dfs_toposort(v, visited);
            }
        }
        prob.tasks_toposorted.push_back(u);
    }

    void calculate_toposort() {
        bool visited[maxn];
        memset(visited, 0, sizeof(visited));
        for (int i = 1; i <= prob.N; ++i) {
            if (!visited[i] && cando[i]) {
                dfs_toposort(i, visited);
            }
        }
        reverse(prob.tasks_toposorted.begin(), prob.tasks_toposorted.end());
        for (int i = 0; i < prob.tasks_toposorted.size(); ++i) {
            task_id_in_toposort_list[prob.tasks_toposorted[i]] = i;
        }
    }

  public:
    Problem prob;
    int cost;
    int penalty;

    int R;                           // total tasks can be done
    vector<int> tasks_of_team[maxm]; // team j do tasks in
                                     // done_by_team[j] sequentially
    bool cando[maxn];                // cando[i] = true if task i can be done
    int start_time[maxn]; // start_time[i] = time to start doing task i
    int task_id_in_toposort_list[maxn]; // task_id_in_toposort_list[i] =
                                        // position of task i in
                                        // tasks_toposorted
    int id_done_by_team[maxn]; // id_done_by_team[i] = order of task i done by
                               // the team
    int team_do[maxn];         // team_do[i] = team do task i

    Solution(Problem prob) {
        this->cost = VERY_BIG_NUMBER;
        this->penalty = -VERY_BIG_NUMBER;
        this->prob = prob;

        for (int i = 1; i <= prob.M; ++i) {
            this->tasks_of_team[i] = vector<int>(0);
        }
        memset(this->tasks_of_team, 0, sizeof(tasks_of_team));
        memset(this->cando, 0, sizeof(cando));
    }

    void update_id_done_by_team() {
        for (int i = 1; i <= prob.M; ++i) {
            for (int j = 0; j < tasks_of_team[i].size(); ++j) {
                int task = tasks_of_team[i][j];
                id_done_by_team[task] = j;
                team_do[task] = i;
            }
        }
    }

    void init() {
        update_cando();
        update_must_before();
        calculate_toposort();

        // assign all task for team 1
        tasks_of_team[1] = prob.tasks_toposorted;

        recalculate_cost();
        recalculate_penalty();
    }

    void calculate_time() {
        // calculate earliest time to do every tasks and store in time[i]
        // ensure inter-team precedence constraints, ie for every u done by i, v
        // done by j, u must be done before v, i!=j,
        // then time[u] + d[u] <= time[v]
        // ASSUME we have a toposorted list
        update_id_done_by_team();
        bool visited[maxn];
        memset(visited, 0, sizeof(visited));
        memset(start_time, 0, sizeof(start_time));
        priority_queue<pair<int, int>, vector<pair<int, int>>,
                       greater<pair<int, int>>>
            pq; // {id_in_toposorted_list, id}
        for (int i = 1; i <= prob.M; ++i) {
            if (tasks_of_team[i].size() > 0) {
                int task_id = tasks_of_team[i][0];
                pq.push({task_id_in_toposort_list[task_id], task_id});
                start_time[task_id] = prob.s[i];
            }
        }
        while (!pq.empty()) {
            auto x = pq.top();
            pq.pop();
            int task_id = x.second;
            int idd_done_by_team = id_done_by_team[task_id];
            int team = team_do[task_id];
            if (visited[task_id])
                continue;
            else
                visited[task_id] = true;
            if (idd_done_by_team < tasks_of_team[team].size() - 1) {
                int next_task_by_same_team =
                    tasks_of_team[team][idd_done_by_team + 1];
                start_time[next_task_by_same_team] =
                    max(start_time[next_task_by_same_team],
                        start_time[task_id] + prob.d[task_id]);
                pq.push({task_id_in_toposort_list[next_task_by_same_team],
                         next_task_by_same_team});
            }
            for (int i = 0; i < prob.succ[task_id].size(); ++i) {
                int v = prob.succ[task_id][i];
                if (cando[v] && team_do[task_id] != team_do[v]) {
                    // ignore the tasks done by same team
                    start_time[v] = max(start_time[v],
                                        start_time[task_id] + prob.d[task_id]);
                    pq.push({task_id_in_toposort_list[v], v});
                }
            }
        }
    }

    bool check_feasible() {
        // for every team:
        // check team i can do tasks in team[i]
        // done_by_team[i] is sequentially done, satisfy precedence constraints

        // not checking inter-team precedence constraints
        for (int i = 1; i <= prob.M; ++i) {
            // check all pairs of tasks done by team i
            for (int j = 0; j < tasks_of_team[i].size(); ++j) {
                int u = tasks_of_team[i][j];
                for (int k = j + 1; k < tasks_of_team[i].size(); ++k) {
                    int v = tasks_of_team[i][k];
                    if (prob.get_must_before(u, v)) {
                        return false;
                    }
                }
            }
            // check if if team i can do every task in done_by_team[i]
            for (int j = 0; j < tasks_of_team[i].size(); ++i) {
                int u = tasks_of_team[i][j];
                if (prob.c[u][i] == 0) {
                    return false;
                }
            }
        }

        return true;
    }

    void recalculate_cost() {
        // earliest time to finish
        calculate_time();
        // cost = max(time[i] + d[i]) for every i
        cost = 0;
        for (int i = 1; i <= prob.N; ++i) {
            cost = max(cost, start_time[i] + prob.d[i]);
        }
    }

    void recalculate_penalty() {
        // number of faults * PENALTY_PER_FAULT
        // fault: not include inter-team precedence constraints (ie, u done by
        // i, v done by j, if u must done before v then start_time[u] + d[u] <=
        // start_time[v])
        int n_fault = 0;
        for (int i = 1; i <= prob.M; ++i) {
            // check all pairs of tasks done by team i
            for (int j = 0; j < tasks_of_team[i].size(); ++j) {
                int u = tasks_of_team[i][j];
                for (int k = j + 1; k < tasks_of_team[i].size(); ++k) {
                    int v = tasks_of_team[i][k];
                    if (prob.get_must_before(u, v)) {
                        ++n_fault;
                    }
                }
            }
            // check if if team i can do every task in done_by_team[i]
            for (int j = 0; j < tasks_of_team[i].size(); ++i) {
                int u = tasks_of_team[i][j];
                if (prob.c[u][i] == 0) {
                    ++n_fault;
                }
            }
        }
        penalty = n_fault * PENALTY_PER_FAULT;
    }

    void print_debug() {
#ifndef NDEBUG
        cout << "cost: " << cost << endl;
        cout << "penalty: " << penalty << endl;
#endif
    }

    void print_answer() {
        cout << R << endl;
        for (int i = 1; i <= prob.N; ++i) {
            cout << i << " " << team_do[i] << " " << start_time[i] << endl;
        }
    }
};

class NeighborOperator {
  public:
    virtual void update_cost(Solution &sol) {}
    virtual void update_penalty(Solution &sol) {}
    virtual Solution find_best_neighbor(Solution sol) { return sol; }
};

class InterSwapOperator : public NeighborOperator {
  public:
    InterSwapOperator() {}

    int most_violate(Solution &sol, int team) {
        int max_violate = -1;
        vector<int> most_violate_pos;
        for (int i = 0; i < sol.tasks_of_team[team].size(); ++i) {
            int n_violate = 0;
            int u = sol.tasks_of_team[team][i];
            if (prob.c[i][team] == 0)
                n_violate++;
            for (int j = 0; j < sol.tasks_of_team[team].size(); ++j) {
                int v = sol.tasks_of_team[team][j];
                if (sol.prob.get_must_before(u, v))
                    n_violate++;
            }
            if (n_violate > max_violate) {
                max_violate = n_violate;
                most_violate_pos.clear();
                most_violate_pos.push_back(i);
            } else if (n_violate == max_violate) {
                most_violate_pos.push_back(i);
            }
        }
        return most_violate_pos[rand() % (most_violate_pos.size())];
    }

    int least_violate_insert_position(Solution &sol, int task, int team) {
        int n_violate_at_0 = 0; // n_violate if insert at position 0
        for (int i = 0; i < sol.tasks_of_team[team].size(); ++i) {
            int u = sol.tasks_of_team[team][i];
            if (sol.prob.get_must_before(u, task))
                n_violate_at_0++;
        }
        int least_violate = n_violate_at_0;
        vector<int> least_violate_pos;
        least_violate_pos.push_back(0);
        for (int i = 0; i < sol.tasks_of_team[team].size(); ++i) {
            int u = sol.tasks_of_team[team][i];
            if (sol.prob.get_must_before(u, task)) {
                n_violate_at_0--;
                if (n_violate_at_0 < least_violate) {
                    least_violate = n_violate_at_0;
                    least_violate_pos.clear();
                    least_violate_pos.push_back(i + 1);
                } else if (n_violate_at_0 == least_violate) {
                    least_violate_pos.push_back(i + 1);
                }
            } else if (sol.prob.get_must_before(task, u)) {
                n_violate_at_0++;
            }
        }
        return least_violate_pos[rand() % (least_violate_pos.size())];
    }

    Solution update_sol(Solution sol, int from_team, int to_team) {
        Solution newsol = sol;
        int which_to_take = most_violate(sol, from_team);
        int task_to_take = newsol.tasks_of_team[from_team][which_to_take];
        int where_to_insert =
            least_violate_insert_position(sol, task_to_take, to_team);
        // int where_to_insert = 0;
        newsol.tasks_of_team[from_team].erase(
            newsol.tasks_of_team[from_team].begin() + which_to_take);
        newsol.tasks_of_team[to_team].insert(
            newsol.tasks_of_team[to_team].begin() + where_to_insert,
            task_to_take);
        return newsol;
    }

    virtual void update_cost(Solution &sol) { sol.recalculate_cost(); }

    virtual void update_penalty(Solution &sol) { sol.recalculate_penalty(); }

    virtual Solution find_best_neighbor(Solution sol) {
        cout << "inter\n";
        int best_cost_penalty = VERY_BIG_NUMBER;
        vector<Solution> best_sols;
        for (int i = 1; i <= prob.M; ++i) {
            if (sol.tasks_of_team[i].size() == 0)
                continue;
            for (int j = 1; j <= prob.M; ++j) {
                if (j == i)
                    continue;
                Solution newsol = update_sol(sol, i, j);
                update_cost(newsol);
                update_penalty(newsol);
                int cost_penalty = newsol.cost + newsol.penalty;
                if (cost_penalty < best_cost_penalty) {
                    best_cost_penalty = cost_penalty;
                    best_sols.clear();
                    best_sols.push_back(newsol);
                } else if (cost_penalty == best_cost_penalty) {
                    best_sols.push_back(newsol);
                }
            }
        }
        return best_sols[rand() % (best_sols.size())];
    }
};

class IntraSwapOperator : public NeighborOperator {
  public:
    IntraSwapOperator() {}

    int most_violate(Solution &sol, int team) {
        int max_violate = -1;
        vector<int> most_violate_pos;
        for (int i = 0; i < sol.tasks_of_team[team].size(); ++i) {
            int n_violate = 0;
            int u = sol.tasks_of_team[team][i];
            if (prob.c[i][team] == 0)
                n_violate++;
            for (int j = 0; j < sol.tasks_of_team[team].size(); ++j) {
                int v = sol.tasks_of_team[team][j];
                if (sol.prob.get_must_before(u, v))
                    n_violate++;
            }
            if (n_violate > max_violate) {
                max_violate = n_violate;
                most_violate_pos.clear();
                most_violate_pos.push_back(i);
            } else if (n_violate == max_violate) {
                most_violate_pos.push_back(i);
            }
        }
        return most_violate_pos[rand() % (most_violate_pos.size())];
    }

    int least_violate_insert_position(Solution &sol, int task, int team) {
        int n_violate_at_0 = 0; // n_violate if insert at position 0
        for (int i = 0; i < sol.tasks_of_team[team].size(); ++i) {
            int u = sol.tasks_of_team[team][i];
            if (sol.prob.get_must_before(u, task))
                n_violate_at_0++;
        }
        int least_violate = n_violate_at_0;
        vector<int> least_violate_pos;
        least_violate_pos.push_back(0);
        for (int i = 0; i < sol.tasks_of_team[team].size(); ++i) {
            int u = sol.tasks_of_team[team][i];
            if (sol.prob.get_must_before(u, task)) {
                n_violate_at_0--;
                if (n_violate_at_0 < least_violate) {
                    least_violate = n_violate_at_0;
                    least_violate_pos.clear();
                    least_violate_pos.push_back(i + 1);
                } else if (n_violate_at_0 == least_violate) {
                    least_violate_pos.push_back(i + 1);
                }
            } else if (sol.prob.get_must_before(task, u)) {
                n_violate_at_0++;
            }
        }
        return least_violate_pos[rand() % (least_violate_pos.size())];
    }


    Solution update_sol(Solution sol, int team, int id_to_move,
                        int id_to_come) {
        Solution newsol = sol;
        int task_to_move = newsol.tasks_of_team[team][id_to_move];
        newsol.tasks_of_team[team].erase(newsol.tasks_of_team[team].begin() +
                                         id_to_move);
        newsol.tasks_of_team[team].insert(
            newsol.tasks_of_team[team].begin() + id_to_come, task_to_move);
        return newsol;
    }

    virtual void update_cost(Solution &sol) { sol.recalculate_cost(); }

    virtual void update_penalty(Solution &sol) { sol.recalculate_penalty(); }

    virtual Solution find_best_neighbor(Solution sol) {
        cout << "intra\n";
        for (int team = 1; team <= prob.M; ++team) {
            int best_cost_penalty = VERY_BIG_NUMBER;
            vector<Solution> best_sols;
            best_sols.push_back(sol);
            for (int i = 0; i < sol.tasks_of_team[team].size(); ++i) {
                for (int j = 0; j < sol.tasks_of_team[team].size(); ++j) {
                    if (j == i)
                        continue;
                    Solution newsol = update_sol(sol, team, i, j);
                    update_cost(newsol);
                    update_penalty(newsol);
                    int cost_penalty = newsol.cost + newsol.penalty;
                    if (cost_penalty < best_cost_penalty) {
                        best_cost_penalty = cost_penalty;
                        best_sols.clear();
                        best_sols.push_back(newsol);
                    } else if (cost_penalty == best_cost_penalty) {
                        best_sols.push_back(newsol);
                    }
                }
            }
            sol = best_sols[rand() % (best_sols.size())];
        }
        return sol;
    }
};

class LocalSearch {
  public:
    vector<NeighborOperator *> neighborOperators;
    vector<double> cmf;

    LocalSearch() {}

    void addNeighborOperator(NeighborOperator &neighborOperator, double prob) {
        assert(prob > 0 && prob < 1);
        neighborOperators.push_back(&neighborOperator);
        for (int i = 0; i < cmf.size(); ++i) {
            cmf[i] *= (1 - prob);
            cmf[i] += 1e-5;
        }
        cmf.push_back(1 + 1e-5);
    }

    Solution search_with_time(Solution sol) {
        int patient = 0;
        Solution best_sol = sol;
        time_t start_time = time(nullptr);
        while (time(nullptr) - start_time < MAXTIME) {
            double r = ((double)rand() / (RAND_MAX));
            int id = upper_bound(cmf.begin(), cmf.end(), r) - cmf.begin();
            Solution neighbor = neighborOperators[id]->find_best_neighbor(sol);
            sol = neighbor;
            // sol.print_answer();
            sol.print_debug();
            if (sol.cost + sol.penalty <= best_sol.cost + best_sol.penalty) {
                best_sol = sol;
                patient = 0;
            } else {
                patient++;
                if (patient >= MAX_PATIENT && best_sol.penalty == 0)
                    break;
            }
        }

        return best_sol;
    }
};

int main() {
    ios::sync_with_stdio(false);
    cin.tie(0);
    srand(time(NULL));

    prob.input();
    Solution sol(prob);
    sol.init();
    InterSwapOperator interSwapOperator;
    IntraSwapOperator intraSwapOperator;
    LocalSearch localsearch;
    localsearch.addNeighborOperator(interSwapOperator, 0.8);
    localsearch.addNeighborOperator(intraSwapOperator, 0.2);
    sol = localsearch.search_with_time(sol);
    sol.print_answer();
    sol.print_debug();

    return 0;
}
