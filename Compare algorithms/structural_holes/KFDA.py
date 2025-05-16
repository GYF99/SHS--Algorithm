import numpy as np
import networkx as nx

def community_detection(G):
    # 使用社区检测算法将网络划分为社区
    return nx.community.greedy_modularity_communities(G)

def network_embedding(G, communities, walk_length=10):
    # Using a community-aware embedding method to embed each node into a low-dimensional vector space

    def random_walk(G, start_node, length):
        visited = set()
        node_list = [start_node]
        for _ in range(length):
            node = node_list[-1]
            neighbors = list(G.neighbors(node))
            if neighbors:
                neighbor = np.random.choice(neighbors)  # Randomly choose a neighbor node
                visited.add(neighbor)
                node_list.append(neighbor)
            else:
                break  # If there are no neighbor nodes, end the random walk
        return node_list

    # Find the maximum community size for padding or truncation
    max_community_size = max(len(community) for community in communities)

    nodes = []
    for community in communities:
        for node in community:
            node_list = random_walk(G, node, walk_length)
            # Convert node names to numeric indices
            node_list_numeric = [list(G.nodes).index(n) for n in node_list]
            # Padding or truncating the node sequence to have the same length
            node_list_numeric += [0] * (walk_length - len(node_list_numeric))
            nodes.append(node_list_numeric)

    X = np.array(nodes)
    return X


def harmonic_centrality(X, Q):
    # Print shapes for debugging
    print("X shape:", X.shape)
    print("Q shape:", Q.shape)

    # 计算每个节点的谐振中心性得分
    return np.sum(X * Q * X, axis=1) / np.sum(Q, axis=1)



def key_figure_detection(G, communities, X, Q):
    # 根据节点的谐振中心性得分对节点进行排名，并应用阈值来识别顶级节点为关键人物
    scores = harmonic_centrality(X, Q)
    return np.argsort(scores)[::-1][:int(0.1 * G.number_of_nodes())]


def harmonic_centrality(X, Q):
    # 计算每个节点的谐振中心性得分
    # Assuming Q is the modularity, you may need to adjust the formula based on your requirements
    return np.sum(X * Q * X, axis=1) / (np.sum(Q) + 1e-8)  # Adding a small epsilon to avoid division by zero


def main():
    # 读取网络数据
    G = nx.read_edgelist("GRQC.txt", create_using=nx.Graph, delimiter=" ")

    # 进行社区检测
    communities = community_detection(G)

    # 进行网络嵌入
    X = network_embedding(G, communities)

    # 计算谐振中心性得分
    Q = nx.community.modularity(G, communities)

    # 识别关键人物
    key_figures = key_figure_detection(G, communities, X, Q)

    # 打印关键人物
    print(key_figures)


if __name__ == "__main__":
    main()








