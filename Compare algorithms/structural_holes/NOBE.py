import easygraph as eg
import numpy as np
import networkx as nx
import igraph as ig
from easygraph.utils import *


__all__ = ["NOBE_SH", "NOBE_GA_SH"]


@not_implemented_for("multigraph")
def NOBE_SH(G, K, topk):
    """detect SH spanners via NOBE[1].

    Parameters
    ----------
    G : easygraph.Graph
        An unweighted and undirected graph.

    K : int
        Embedding dimension k

    topk : int
        top - k structural hole spanners

    Returns
    -------
    SHS : list
        The top-k structural hole spanners.

    Examples
    --------
    >>> NOBE_SH(G,K=8,topk=5)

    References
    ----------
    .. [1] https://www.researchgate.net/publication/325004496_On_Spectral_Graph_Embedding_A_Non-Backtracking_Perspective_and_Graph_Approximation

    """
    from sklearn.cluster import KMeans

    Y = eg.graph_embedding.NOBE(G, K)
    dict = {}
    a = 0
    for i in G.nodes:
        dict[i] = a
        a += 1
    if isinstance(Y[0, 0], complex):
        Y = abs(Y)
    kmeans = KMeans(n_clusters=K, random_state=0).fit(Y)
    com = {}
    cluster = {}
    for i in dict:
        com[i] = kmeans.labels_[dict[i]]
    for i in com:
        if com[i] in cluster:
            cluster[com[i]].append(i)
        else:
            cluster[com[i]] = []
            cluster[com[i]].append(i)
    vector = {}
    for i in dict:
        vector[i] = Y[dict[i]]
    rds = RDS(com, cluster, vector, K)
    rds_sort = sorted(rds.items(), key=lambda d: d[1], reverse=True)
    SHS = list()
    a = 0
    for i in rds_sort:
        SHS.append(i[0])
        a += 1
        if a == topk:
            break
    return SHS


@not_implemented_for("multigraph")
def NOBE_GA_SH(G, K, topk):
    """detect SH spanners via NOBE-GA[1].

    Parameters
    ----------
    G : easygraph.Graph
        An unweighted and undirected graph.

    K : int
        Embedding dimension k

    topk : int
        top - k structural hole spanners

    Returns
    -------
    SHS : list
        The top-k structural hole spanners.

    Examples
    --------
    >>> NOBE_GA_SH(G,K=8,topk=5)

    References
    ----------
    .. [1] https://www.researchgate.net/publication/325004496_On_Spectral_Graph_Embedding_A_Non-Backtracking_Perspective_and_Graph_Approximation

    """
    from sklearn.cluster import KMeans

    Y = eg.NOBE_GA(G, K)
    if isinstance(Y[0, 0], complex):
        Y = abs(Y)
    kmeans = KMeans(n_clusters=K, random_state=0).fit(Y)
    com = {}
    cluster = {}
    a = 0
    for i in G.nodes:
        com[i] = kmeans.labels_[a]
        a += 1
    for i in com:
        if com[i] in cluster:
            cluster[com[i]].append(i)
        else:
            cluster[com[i]] = []
            cluster[com[i]].append(i)
    vector = {}
    a = 0
    for i in G.nodes:
        vector[i] = Y[a]
        a += 1
    rds = RDS(com, cluster, vector, K)
    rds_sort = sorted(rds.items(), key=lambda d: d[1], reverse=True)
    SHS = list()
    a = 0
    for i in rds_sort:
        SHS.append(i[0])
        a += 1
        if a == topk:
            break
    return SHS


def RDS(com, cluster, vector, K):
    rds = {}
    Uc = {}
    Rc = {}
    for i in cluster:
        sum_vec = np.zeros(K)
        for j in cluster[i]:
            sum_vec += vector[j]
        Uc[i] = sum_vec / len(cluster[i])
    for i in cluster:
        sum_dist = 0
        for j in cluster[i]:
            sum_dist += np.linalg.norm(vector[j] - Uc[i])
        Rc[i] = sum_dist
    for i in com:
        maxx = 0
        fenzi = np.linalg.norm(vector[i] - Uc[com[i]]) / Rc[com[i]]
        for j in cluster:
            fenmu = np.linalg.norm(vector[i] - Uc[j]) / Rc[j]
            if maxx < fenzi / fenmu:
                maxx = fenzi / fenmu
        rds[i] = maxx
    return rds


if __name__ == "__main__":
    G = eg.datasets.get_graph_karateclub()
    # ## 真实网络
    # path = r"./"
    # ## 真实网络
    # karate_network = path + r'karate.txt'  #
    # zhang_network = path + r'zhang.txt'  #
    # football_network = path + r'football.txt'  #
    # dolphins_network = path + r'dolphins.txt'  #
    # # 选择网络
    # # real
    # # network = dolphins_network
    # # network = football_network
    # # network = zhang_network
    # network = karate_network
    # network_name = 'karate'
    # # network_name = 'zhang'
    # # network_name = 'football'
    # # network_name = 'dolphins'
    # groundtruth_path = path + "/real/" + network_name + '_groundtruth.txt'
    # with open(groundtruth_path, mode='r', encoding='UTF-8') as f:
    #     real_mem = list(map(int, f.read().splitlines()))
    # # 获取网络数据中的边列表，并根据其使用igraph创建网络
    # G1 = nx.read_edgelist(network, create_using=nx.Graph())
    # # 转换 G1 为无向图（在 NetworkX 3.1 中已默认读取无向图）
    # # 所以这一步在 NetworkX 3.1 中其实不是必需的
    # # G1 = G1.to_undirected()
    # # 从 G1 转换为 igraph 中的图对象 Gi
    # Gi = ig.Graph.Read_Edgelist(network)
    # # 由于 igraph 中的节点必须是整数，所以需要将 G1 中的节点转换为整数
    # # 并在 Gi 中保留相应的子图
    # Gi = Gi.subgraph(map(int, G1.nodes()))
    # # 将 igraph 中的图对象 Gi 转换回 NetworkX 图对象 G
    # # 并转换为无向图（在 NetworkX 3.1 中默认读取无向图）
    #
    # G = nx.Graph(Gi.get_edgelist())
    # G.cflag = 0
    # mapping = {node: i for i, node in enumerate(G.nodes())}
    # G = nx.relabel_nodes(G, mapping)

    print(NOBE_SH(G, K=2, topk=3))
    # print(NOBE_GA_SH(G, K=2, topk=3))
