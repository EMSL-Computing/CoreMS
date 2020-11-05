import numpy as np
import scipy.stats
'''exploratory module based on  Yuanyue Li code at 
TODO add GitHub and Paper here'''

def entropy_distance(p, q):
    merged = p + q
    entropy_increase = 2 * scipy.stats.entropy(merged) - scipy.stats.entropy(p) - scipy.stats.entropy(q)
    return entropy_increase

def _weight_intensity_for_entropy(x):
    if sum(x) > 0:
        WEIGHT_START = 0.25
        WEIGHT_SLOPE = 0.5

        entropy_x = scipy.stats.entropy(x)
        weight = WEIGHT_START + WEIGHT_SLOPE * entropy_x
        x = np.power(x, weight)
        x = x / sum(x)
        return x


def weighted_entropy_distance(p, q):
    p = _weight_intensity_for_entropy(p)
    q = _weight_intensity_for_entropy(q)

    merged = p + q
    entropy_increase = 2 * scipy.stats.entropy(merged) - scipy.stats.entropy(p) - scipy.stats.entropy(q)
    return entropy_increase


def chebyshev_distance(p, q):
    r"""
    Chebyshev distance:

    .. math::

        \underset{i}{\max}{(|P_{i}\ -\ Q_{i}|)}
    """
    return np.max(np.abs(p - q))


def squared_euclidean_distance(p, q):
    r"""
    Squared Euclidean distance:

    .. math::

        \sum(P_{i}-Q_{i})^2
    """
    return np.sum(np.power(p - q, 2))


def fidelity_distance(p, q):
    r"""
    Fidelity distance:

    .. math::

        1-\sum\sqrt{P_{i}Q_{i}}
    """
    return 1 - np.sum(np.sqrt(p * q))


def matusita_distance(p, q):
    r"""
    Matusita distance:

    .. math::

        \sqrt{\sum(\sqrt{P_{i}}-\sqrt{Q_{i}})^2}
    """
    return np.sqrt(np.sum(np.power(np.sqrt(p) - np.sqrt(q), 2)))


def squared_chord_distance(p, q):
    r"""
    Squared-chord distance:

    .. math::

        \sum(\sqrt{P_{i}}-\sqrt{Q_{i}})^2
    """
    return np.sum(np.power(np.sqrt(p) - np.sqrt(q), 2))


def bhattacharya_1_distance(p, q):
    r"""
    Bhattacharya 1 distance:

    .. math::

        (\arccos{(\sum\sqrt{P_{i}Q_{i}})})^2
    """
    s = np.sum(np.sqrt(p * q))
    # TODO:Fix this!
    if s > 1:
        if s > 1 + 1e-6:
            print("Error in calculating Bhattacharya 1 distance, got arccos {}".format(s))
        s = 1
    return np.power(np.arccos(s), 2)


def bhattacharya_2_distance(p, q):
    r"""
    Bhattacharya 2 distance:

    .. math::

        -\ln{(\sum\sqrt{P_{i}Q_{i}})}
    """
    s = np.sum(np.sqrt(p * q))
    if s == 0:
        return np.inf
    else:
        return -np.log(s)


def harmonic_mean_distance(p, q):
    r"""
    Harmonic mean distance:

    .. math::

        1-2\sum(\frac{P_{i}Q_{i}}{P_{i}+Q_{i}})
    """
    return 1 - 2 * np.sum(p * q / (p + q))


def pearson_chi_squared_distance(p, q):
    r"""
    Pearson χ2 distance:

    .. math::

        \sum\frac{(P_{i}-Q_{i})^2}{Q_{i}}
    """
    return np.sum(np.power(p - q, 2) / q)


def neyman_chi_squared_distance(p, q):
    r"""
    Neyman χ2 distance:

    .. math::

        \sum\frac{(P_{i}-Q_{i})^2}{P_{i}}
    """
    return np.sum(np.power(p - q, 2) / p)


def probabilistic_symmetric_chi_squared_distance(p, q):
    r"""
    Probabilistic symmetric χ2 distance:

    .. math::

        \frac{1}{2} \times \sum\frac{(P_{i}-Q_{i}\ )^2}{P_{i}+Q_{i}\ }
    """
    return 1 / 2 * np.sum(np.power(p - q, 2) / (p + q))


def topsoe_distance(p, q):
    r"""
    Topsøe distance:

    .. math::

        \sum{(P_{i}ln\frac{P_{i}}{Z_i}+Q_{i}ln\frac{Q_{i}}{Z_i}),\ \ \ Z_i=\frac{1}{2}(P_{i}+Q_{i})}
    """
    z = 1 / 2 * (p + q)
    z[z == 0] = 1
    pz = p / z
    qz = q / z
    pz[p == 0] = 1
    qz[q == 0] = 1
    return np.sum(p * np.log(pz) + q * np.log(qz))


def chernoff_distance(p, q):
    r"""
    Chernoff distance:

    .. math::

        \max{(-ln\sum(P_{i}^tQ_{i}^{1-t})^{1-t})},\ t=0.1,\ 0\le\ t<1
    """
    t = 0.1
    return np.max(-np.log(
        np.sum(np.power(np.power(p, t) * np.power(q, 1 - t), 1 - t))))


def ruzicka_distance(p, q):
    r"""
    Ruzicka distance:

    .. math::

        \frac{\sum{|P_{i}-Q_{i}|}}{\sum{\max(P_{i},Q_{i})}}
    """
    dist = np.sum(np.abs(p - q)) / np.sum(np.maximum(p, q))
    return dist


def roberts_distance(p, q):
    r"""
    Roberts distance:

    .. math::

        1-\sum\frac{(P_{i}+Q_{i})\frac{\min{(P_{i},Q_{i})}}{\max{(P_{i},Q_{i})}}}{\sum(P_{i}+Q_{i})}
    """
    return 1 - np.sum((p + q) * np.minimum(p, q) / np.maximum(p, q) / np.sum(p + q))


def intersection_distance(p, q):
    r"""
    Intersection distance:

    .. math::

        1-\frac{\sum\min{(P_{i},Q_{i})}}{\min(\sum{P_{i},\sum{Q_{i})}}}
    """
    return 1 - np.sum(np.minimum(p, q)) / min(np.sum(p), np.sum(q))


def motyka_distance(p, q):
    r"""
    Motyka distance:

    .. math::

        -\frac{\sum\min{(P_{i},Q_{i})}}{\sum(P_{i}+Q_{i})}
    """
    dist = np.sum(np.minimum(p, q)) / np.sum(p + q)
    return -dist


def canberra_distance(p, q):
    r"""
    Canberra distance:

    .. math::

        \sum\frac{|P_{i}-Q_{i}|}{|P_{i}|+|Q_{i}|}
    """
    return np.sum(np.abs(p - q) / (np.abs(p) + np.abs(q)))


def kulczynski_1_distance(p, q):
    r"""
    Kulczynski 1 distance:

    .. math::

        \frac{\sum{|P_i}-Q_i|}{\sum m\ i\ n\ (P_i,Q_i)}
    """
    return np.sum(np.abs(p - q)) / np.sum(np.minimum(p, q))


def baroni_urbani_buser_distance(p, q):
    r"""
    Baroni-Urbani-Buser distance:

    .. math::

        1-\frac{\sum\min{(P_i,Q_i)}+\sqrt{\sum\min{(P_i,Q_i)}\sum(\max{(P)}-\max{(P_i,Q_i)})}}{\sum{\max{(P_i,Q_i)}+\sqrt{\sum{\min{(P_i,Q_i)}\sum(\max{(P)}-\max{(P_i,Q_i)})}}}}
    """
    if np.max(p) < np.max(q):
        p, q = q, p
    d1 = np.sqrt(np.sum(np.minimum(p, q) * np.sum(max(p) - np.maximum(p, q))))
    return 1 - (np.sum(np.minimum(p, q)) + d1) / (np.sum(np.maximum(p, q)) + d1)


def penrose_size_distance(p, q):
    r"""
    Penrose size distance:

    .. math::

        \sqrt N\sum{|P_i-Q_i|}
    """
    n = np.sum(p > 0)
    return np.sqrt(n) * np.sum(np.abs(p - q))


def mean_character_distance(p, q):
    r"""
    Mean character distance:

    .. math::

        \frac{1}{N}\sum{|P_i-Q_i|}
    """
    n = np.sum(p > 0)
    return 1 / n * np.sum(np.abs(p - q))


def lorentzian_distance(p, q):
    r"""
    Lorentzian distance:

    .. math::

        \sum{\ln(1+|P_i-Q_i|)}
    """
    return np.sum(np.log(1 + np.abs(p - q)))


def penrose_shape_distance(p, q):
    r"""
    Penrose shape distance:

    .. math::

        \sqrt{\sum((P_i-\bar{P})-(Q_i-\bar{Q}))^2}
    """
    p_avg = np.mean(p)
    q_avg = np.mean(q)
    return np.sqrt(np.sum(np.power((p - p_avg) - (q - q_avg), 2)))


def clark_distance(p, q):
    r"""
    Clark distance:

    .. math::

        (\frac{1}{N}\sum(\frac{P_i-Q_i}{|P_i|+|Q_i|})^2)^\frac{1}{2}
    """
    n = np.sum(p > 0)
    return np.sqrt(1 / n * np.sum(np.power((p - q) / (np.abs(p) + np.abs(q)), 2)))


def hellinger_distance(p, q):
    r"""
    Hellinger distance:

    .. math::

        \sqrt{2\sum(\sqrt{\frac{P_i}{\bar{P}}}-\sqrt{\frac{Q_i}{\bar{Q}}})^2}
    """
    p_avg = np.mean(p)
    q_avg = np.mean(q)
    return np.sqrt(2 * np.sum(np.power(np.sqrt(p / p_avg) - np.sqrt(q / q_avg), 2)))


def whittaker_index_of_association_distance(p, q):
    r"""
    Whittaker index of association distance:

    .. math::

        \frac{1}{2}\sum|\frac{P_i}{\bar{P}}-\frac{Q_i}{\bar{Q}}|
    """
    p_avg = np.mean(p)
    q_avg = np.mean(q)
    return 1 / 2 * np.sum(np.abs(p / p_avg - q / q_avg))


def symmetric_chi_squared_distance(p, q):
    r"""
    Symmetric χ2 distance:

    .. math::

        \sqrt{\sum{\frac{\bar{P}+\bar{Q}}{N(\bar{P}+\bar{Q})^2}\frac{(P_i\bar{Q}-Q_i\bar{P})^2}{P_i+Q_i}\ }}
    """
    p_avg = np.mean(p)
    q_avg = np.mean(q)
    n = np.sum(p > 0)

    d1 = (p_avg + q_avg) / (n * np.power(p_avg + q_avg, 2))
    return np.sqrt(d1 * np.sum(np.power(p * q_avg - q * p_avg, 2) / (p + q)))

def similarity_index_distance(p, q):
    r"""
    Similarity Index Distance:

    .. math::

        \sqrt{\frac{\sum\{\frac{P_i-Q_i}{Q_i}\}^2}{N}}
    """
    n = np.sum(p > 0)
    return np.sqrt(1 / n * np.sum(np.power((p - q) / q, 2)))


def improved_similarity_distance(p, q):
    r"""
    Improved Similarity Index:

    .. math::

        \sqrt{\frac{1}{N}\sum\{\frac{P_i-Q_i}{P_i+Q_i}\}^2}
    """
    n = np.sum(p > 0)
    return np.sqrt(1 / n * np.sum(np.power((p - q) / (p + q), 2)))


def absolute_value_distance(p, q):
    r"""
    Absolute Value Distance:

    .. math::

        \frac { \sum(|Q_i-P_i|)}{\sum P_i}

    """
    dist = np.sum(np.abs(q - p)) / np.sum(p)
    return dist

def spectral_contrast_angle_distance(p, q):
    r"""
    Spectral Contrast Angle:

    .. math::

        1 - \frac{\sum{Q_iP_i}}{\sqrt{\sum Q_i^2\sum P_i^2}}
    """
    return 1 - np.sum(q * p) / \
           np.sqrt(np.sum(np.power(q, 2)) * np.sum(np.power(p, 2)))


def wave_hedges_distance(p, q):
    r"""
    Wave Hedges distance:

    .. math::

        \sum\frac{|P_i-Q_i|}{\max{(P_i,Q_i)}}
    """
    return np.sum(np.abs(p - q) / np.maximum(p, q))

def dice_distance(p, q):
    r"""
    Dice distance:

    .. math::

        \frac{\sum(P_i-Q_i)^2}{\sum P_i^2+\sum Q_i^2}
    """
    return np.sum(np.power(p - q, 2)) / \
           (np.sum(np.power(p, 2)) + np.sum(np.power(q, 2)))


def inner_product_distance(p, q):
    r"""
    Inner Product distance:

    .. math::

        1-\sum{P_iQ_i}
    """
    return 1 - np.sum(p * q)


def divergence_distance(p, q):
    r"""
    Divergence distance:

    .. math::

        2\sum\frac{(P_i-Q_i)^2}{(P_i+Q_i)^2}
    """
    return 2 * np.sum((np.power(p - q, 2)) / np.power(p + q, 2))


def _chi_squared_distance(p, q):
    r"""
    Additive symmetric χ2 distance:

    .. math::

        \sum\frac{(P_i-Q_i)^2(P_i+Q_i)}{P_iQ_i}
    """
    dist = np.sum(np.power(p - q, 2) * (p + q) / (p * q))
    return dist


def jensen_difference_distance(p, q):
    r"""
    Jensen difference:

    .. math::

        \sum[\frac{1}{2}(P_i\ln{P_i}+Q_i\ln{Q_i})-(\frac{P_i+Q_i}{2})\ln{(\frac{P_i+Q_i}{2})}]
    """
    p_q_avg = (p + q) / 2
    return np.sum(
        1 / 2 * (p * np.log(p) + q * np.log(q)) -
        p_q_avg * np.log(p_q_avg)
    )


def kumar_johnson_distance(p, q):
    r"""
    Kumar-Johnson distance:

    .. math::

        \sum\frac{(P_i^2-Q_i^2)^2}{2(P_iQ_i)^\frac{3}{2}}
    """
    return np.sum(
        np.power(np.power(p, 2) - np.power(q, 2), 2) / \
        (2 * np.power(p * q, 3 / 2))
    )


def avg_l_distance(p, q):
    r"""
    Avg (L1, L∞) distance:

    .. math::

        \frac{1}{2}(\sum|P_i-Q_i|+\underset{i}{\max}{|P_i-Q_i|})
    """
    return 1 / 2 * (np.sum(np.abs(p - q)) + max(np.abs(p - q)))


def vicis_wave_hadges_distance(p, q):
    r"""
    Vicis-Wave Hadges distance:

    .. math::

        \sum\frac{|P_i-Q_i|}{\min{(P_i,\ Q_i)}}
    """
    return np.sum(np.abs(p - q) / np.minimum(p, q))


def vicis_symmetric_chi_squared_1_distance(p, q):
    r"""
    Vicis-Symmetric χ2 1 distance:

    .. math::

        \sum\frac{(P_i-Q_i)^2}{\min{(P_i,Q_i)^2}}
    """
    return np.sum(np.power(p - q, 2) / np.power(np.minimum(p, q), 2))


def vicis_symmetric_chi_squared_2_distance(p, q):
    r"""
    Vicis-Symmetric χ2 2 distance:

    .. math::

        \sum\frac{(P_i-Q_i)^2}{\min{(P_i,Q_i)}}
    """
    return np.sum(np.power(p - q, 2) / np.minimum(p, q))


def vicis_symmetric_chi_squared_3_distance(p, q):
    r"""
    Vicis-Symmetric χ2 3 distance:

    .. math::

        \sum\frac{(P_i-Q_i)^2}{\max{(P_i,Q_i)}}
    """
    return np.sum(np.power(p - q, 2) / np.maximum(p, q))


def max_symmetric_chi_squared_distance(p, q):
    r"""
    Max-Symmetric χ2 distance:

    .. math::

        \max{(\sum\frac{(P_i-Q_i)^2}{P_i},\sum\frac{(P_i-Q_i)^2}{Q_i})}
    """
    return max(np.sum(np.power(p - q, 2) / p), np.sum(np.power(p - q, 2) / q))


def min_symmetric_chi_squared_distance(p, q):
    r"""
    Min-Symmetric χ2 distance:

    .. math::

        \min{(\sum\frac{(P_i-Q_i)^2}{P_i},\sum\frac{(P_i-Q_i)^2}{Q_i})}
    """
    return min(np.sum(np.power(p - q, 2) / p), np.sum(np.power(p - q, 2) / q))



