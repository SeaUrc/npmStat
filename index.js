/* --- DESCRIPTIVE --- */

function min(a, comp = (u, v) => {
    return u < v
}) {
    if (!a.length) {
        throw new Error("Array is empty");
    }
    let mn = a[0];
    for (let i = 0; i < a.length; i++) {
        if (comp(a[i], mn)) {
            mn = a[i];
        }
    }
    return mn
}

function max(a, comp = (u, v) => {
    return u > v
}) {
    if (!a.length) {
        throw new Error("Array is empty");
    }
    let mn = a[0];
    for (let i = 0; i < a.length; i++) {
        if (comp(a[i], mn)) {
            mn = a[i];
        }
    }
    return mn
}

function range(a) {
    if (!a.length) {
        throw new Error("Array is empty");
    }
    return max(a) - min(a);
}

function quartile(a, percentile) {
    if (!a.length) {
        throw new Error("Array is empty");
    }
    if (percentile < 0 || percentile > 100) {
        throw new Error("Percentile out of bounds");
    }

    a.sort((a, b) => a - b);

    let N = a.length;
    let rank = (percentile / 100) * (N); // 0 -> N

    if (rank < 1) {
        return a[0];
    }
    rank--;

    let lowerIndex = Math.floor(rank);
    let upperIndex = Math.ceil(rank);
    let lowerValue = a[lowerIndex];
    let upperValue = a[upperIndex];
    let frac = rank - Math.floor(rank);
    return lowerValue + (upperValue - lowerValue) * frac;
}

function median(a) {
    if (!a.length) {
        throw new Error("Array is empty")
    }
    return quartile(a, 50);
}

function Q1(a) {
    if (!a.length) {
        throw new Error("Array is empty")
    }
    return quartile(a, 25);
}

function Q3(a) {
    if (!a.length) {
        throw new Error("Array is empty")
    }
    return quartile(a, 75);
}

function IQR(a) {
    if (!a.length) {
        throw new Error("Array is empty")
    }
    return Q3(a) - Q1(a);
}

function boxPlot(a, excludeOutliers = false) {
    if (!excludeOutliers) {
        return [min(a), Q1(a), median(a), Q3(a), max(a)]
    }

    let threshold = [Q1(a) - 1.5 * IQR(a), Q3(a) + 1.5 * IQR(a)];
    let b = [];
    let outliers = []
    a.forEach((value) => {
        if (value >= threshold[0] && value <= threshold[1]) {
            b.push(value);
        } else {
            outliers.push(value);
        }
    })
    return [min(b), Q1(b), median(b), Q3(b), max(b), outliers];
}

function sum(a) {
    if (!a.length) {
        throw new Error("Array is empty");
    }
    let sum = 0;
    a.forEach((val) => {
        sum += val;
    })
    return sum;
}

function mean(a) {
    if (!a.length) {
        throw new Error("Array is empty");
    }
    return sum(a) / a.length;
}

function sampleVariance(a) {
    if (!a.length) {
        throw new Error("Array is empty");
    }
    let m = mean(a);
    let n = a.length;
    let rollingSum = 0;
    a.forEach((val) => {
        rollingSum += (val - m) ** 2;
    })
    return rollingSum / (n - 1);
}

function populationVariance(a) {
    if (!a.length) {
        throw new Error("Array is empty");
    }
    let m = mean(a);
    let n = a.length;
    let rollingSum = 0;
    a.forEach((val) => {
        rollingSum += (val - m) ** 2;
    })
    return rollingSum / (n);
}

function sampleStd(a) {
    if (!a.length) {
        throw new Error("Array is empty");
    }
    return Math.sqrt(sampleVariance(a));
}

function populationStd(a) {
    if (!a.length) {
        throw new Error("Array is empty");
    }
    return Math.sqrt(populationVariance(a));
}

function coefficientOfVariation(a) {
    if (!a.length) {
        throw new Error("Array is empty");
    }
    return (populationStd(a) / mean(a));
}

/* --- Probability Distributions --- */


/* --- Correlation and Regression --- */
function pearsonCorrelation(x, y) {
    if (!x.length || !y.length) {
        throw new Error("Empty array");
    }
    if (x.length != y.length) {
        throw new Error("Lengths don't match");
    }
    const mx = mean(x);
    const my = mean(y);
    let num = 0;
    for (let i = 0; i < x.length; i++) {
        num += (x[i] - mx) * (y[i] - my);
    }
    let d1 = 0, d2 = 0;
    for (let i = 0; i < x.length; i++) {
        d1 += (x[i] - mx) ** 2;
        d2 += (y[i] - my) ** 2;
    }
    let d = Math.sqrt(d1 * d2);
    return num / d;
}

/* --- Hypothesis Testing --- */


/* --- RNG --- */

/* --- Matrix Operations --- */

/* --- Misc --- */
function factorial(x) {
    if (!Number.isInteger(x)) {
        throw new Error("Non integer factorial. If this was intentional, use the gamma function")
    }
    if (x == 1) {
        return 1;
    }
    return x * factorial(x - 1);
}

function gamma(z) {
    // Lanczos approximation
    let g = 7;
    let C = [0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];
    if (z < 0.5) {
        return Math.PI / (Math.sin(Math.PI * z) * gamma(1 - z));
    } else {
        z -= 1;
        let x = C[0];
        for (var i = 1; i < g + 2; i++)
            x += C[i] / (z + i);

        let t = z + g + 0.5;
        return Math.sqrt(2 * Math.PI) * Math.pow(t, (z + 0.5)) * Math.exp(-t) * x;
    }
}

function choose(N, K){
    if (!Number.isInteger(N) || !Number.isInteger(K)){
        throw new Error("Non integer in combinatoric");
    }
    return ((factorial(N)) / (factorial(K) * factorial(N-K)));
}

module.exports = {
    min,
    max,
    range,
    quartile,
    median,
    Q1,
    Q3,
    IQR,
    boxPlot,
    sum,
    mean,
    sampleVariance,
    populationVariance,
    sampleStd,
    populationStd,
    coefficientOfVariation,
    pearsonCorrelation,
    factorial,
    gamma,
    choose
}