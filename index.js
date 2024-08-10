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
/*
* Args
* normalpdf(x, mean, standard variation)
* returns the probability density of normal distribution at x
* */
function normalpdf(x, mu, sigma){
    const sqrtTwoPi = Math.sqrt(2 * Math.PI);
    const expInside = -1*((x-mu)**2 / (2*(sigma**2)));
    return (1/(sigma*sqrtTwoPi)) * Math.exp(expInside);
}

/*
* Args
* normalcdf(upper bound, mean, standard variation)
* normalcdf(lower bound, upper bound, mean, standard variation)
* returns the cumulative probability from the lower to upper bound of the normal distribution
* */
function normalcdf(...args){
    let lower, upper, mu, sigma;
    if (args.length === 3) {
        [upper, mu, sigma] = args;
        lower = Number.MIN_VALUE;  // Assume the lower bound is min value
    } else if (args.length === 4) {
        [lower, upper, mu, sigma] = args;
    } else {
        throw new Error("Invalid number of arguments. Expected 3 or 4 arguments.");
    }
    return (1+erf((upper-mu)/(Math.sqrt(2)*sigma)))*0.5  -  (1+erf((lower-mu)/(Math.sqrt(2)*sigma)))*0.5;
}

/*
* Args
* invNorm(probability, mean, standard deviation)
* returns x such that the cumulative probability is equal to the given probability
* */
function invNorm(p, mu=0, sigma=1, tail="lower") {
    if (p < 0 || p > 1) {
        throw new Error("Given probability must be in the range [0, 1]");
    }
    if (p == 0){
        return Number.MIN_VALUE;
    }
    if (p == 1){
        return Number.MAX_VALUE;
    }
    let pPrime = p;
    if (tail.toLowerCase() == "lower"){
        pPrime=p;
    }else if (tail.toLowerCase() == "upper"){
        pPrime = 1-p;
    }else if (tail.toLowerCase() == "center"){
        pPrime = 1 - ((1-p) / 2);
    }else{
        throw new Error("Tail not an available tail. lower, upper, and center are the only given tails");
    }
    let res = mu + sigma * Math.sqrt(2) * invErf(2 * pPrime - 1);
    if (tail.toLowerCase() == "center"){
        return [mu-(res-mu), res];
    }
    return res;
}

/*
* Args
* binompdf(trials, probability of success for each trial, successes)
* returns the probability density of binomial distribution with N, P with K successes
* */
function binompdf(N, P, K){
    if (!Number.isInteger(N) || !Number.isInteger(K)){
        throw new Error("N and K must be integers");
    }
    if (K > N){
        throw new Error("# of successes cannot be higher than number of trials!");
    }
    if (K < 0){
        throw new Error("Cannot have a negative number of successes");
    }
    if (P < 0){
        throw new Error("Cannot have a negative probabilty of success");
    }
    if (P > 1){
        throw new Error("Probability of success cannot be greater than 1");
    }
   return choose(N, K) * (P**K)* ((1-P)**(N-K));
}

/*
* Args
* binomcdf(trials, probability of success for each trial, successes)
* returns the cumulative probabilty of the binomial distribution of K or fewer successes
* */
function binomcdf(N, P, K){
    if (!Number.isInteger(N) || !Number.isInteger(K)){
        throw new Error("N and K must be integers");
    }
    if (K > N){
        throw new Error("# of successes cannot be higher than number of trials!");
    }
    if (K < 0){
        throw new Error("Cannot have a negative number of successes");
    }
    if (P < 0){
        throw new Error("Cannot have a negative probabilty of success");
    }
    if (P > 1){
        throw new Error("Probability of success cannot be greater than 1");
    }
    let cdf = 0;
    for (let i =0; i<=K; i++){
        cdf += binompdf(N, P, i);
    }
    return cdf;
}

/*
* Args
* invBinom(cumulative probability, trials, probability of success in a trial)
* returns the smallest value k such that the cumulative distribution is greater than or equal to the probability given
* */
function invBinom(probability, N, P) {
    if (!Number.isInteger(N)){
        throw new Error("N must be integers");
    }
    if (P < 0){
        throw new Error("Cannot have a negative probabilty of success");
    }
    if (P > 1){
        throw new Error("Probability of success cannot be greater than 1");
    }
    if (P < 0){
        throw new Error("Cannot have a negative cumulative probability");
    }
    if (P > 1){
        throw new Error("Cannot have a cumulative probability greater than 1");
    }
    let cProb = 0;
    for (let K = 0; K <= N; K++) {
        cProb += binompdf(N, P, k);
        if (cProb >= probability) {
            return K;
        }
    }
    // If not found, return n;
    return N;
}

/*
* Args
* tpdf(x, degrees of freedom)
* returns the probability density of the t-distribution with given degrees of freedom at x
* */
function tpdf(x, df) {
    let firstNumerator = gamma((df+1)/2);
    let firstDemon = Math.sqrt(Math.PI*df) * gamma(df/2);
    let second = 1+((x**2)/df);
    return firstNumerator/firstDemon * second**(-1*(df+1)/2);
}

function tcdf(...args){
    let lower, upper, df;
    if (args.length === 2) {
        [upper, df] = args;
        lower = Number.MIN_VALUE;  // Assume the lower bound is min value
    } else if (args.length === 3) {
        [lower, upper, df] = args;
    } else {
        throw new Error("Invalid number of arguments. Expected 2 or 3 arguments.");
    }
    let func = (t) => {return (df/(t**2 + df))};
    let cdfx = (x) => {
        if (x == 0){
            return 0.5;
        }
        if (x < 0){
            return (0.5*regularizedIncompleteBeta(func(x), df/2, 0.5));
        }
        return (1 - 0.5*regularizedIncompleteBeta(func(x), df/2, 0.5))
    };
    return cdfx(upper) - cdfx(lower);
}

function invT(p, df) {
    if (p <= 0.0 || p >= 1.0) {
        throw new Error("p must be between 0 and 1");
    }

    const q = invNorm(p);
    const q2 = q * q;
    const a = (q2 + 1) / (4 * df);
    const b = ((5 * q2 + 16) * q2 + 3) / (96 * df * df);
    const c = (((3 * q2 + 19) * q2 + 17) * q2 - 15) / (384 * df * df * df);
    const d = ((((79 * q2 + 776) * q2 + 1482) * q2 - 1920) * q2 - 945) / (92160 * df * df * df * df);

    const quantile = q * (1 + a + b + c + d);
    return quantile;
}


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

function erf(x) {
    let z;
    const ERF_A = 0.147;
    let the_sign_of_x;
    if (0 == x) {
        the_sign_of_x = 0;
        return 0;
    } else if (x > 0) {
        the_sign_of_x = 1;
    } else {
        the_sign_of_x = -1;
    }

    let one_plus_axsqrd = 1 + ERF_A * x * x;
    let four_ovr_pi_etc = 4 / Math.PI + ERF_A * x * x;
    let ratio = four_ovr_pi_etc / one_plus_axsqrd;
    ratio *= x * -x;
    let expofun = Math.exp(ratio);
    let radical = Math.sqrt(1 - expofun);
    z = radical * the_sign_of_x;
    return z;
}

// https://stackoverflow.com/questions/12556685/is-there-a-javascript-implementation-of-the-inverse-error-function-akin-to-matl
function invErf(x) {
    let z;
    let a = 0.147;
    let the_sign_of_x;
    if (0 == x) {
        the_sign_of_x = 0;
    } else if (x > 0) {
        the_sign_of_x = 1;
    } else {
        the_sign_of_x = -1;
    }

    if (0 != x) {
        let ln_1minus_x_sqrd = Math.log(1 - x * x);
        let ln_1minusxx_by_a = ln_1minus_x_sqrd / a;
        let ln_1minusxx_by_2 = ln_1minus_x_sqrd / 2;
        let ln_etc_by2_plus2 = ln_1minusxx_by_2 + (2 / (Math.PI * a));
        let first_sqrt = Math.sqrt((ln_etc_by2_plus2 * ln_etc_by2_plus2) - ln_1minusxx_by_a);
        let second_sqrt = Math.sqrt(first_sqrt - ln_etc_by2_plus2);
        z = second_sqrt * the_sign_of_x;
    } else { // x is zero
        z = 0;
    }
    return z;
}

function beta(z1, z2){
    return (gamma(z1)*gamma(z2))/gamma(z1+z2);
}

function recurseContinuedFractionIncompleteBeta(m, x, a, b, maxM){
    if (m >= maxM){
        return 1;
    }
    if (m == 0){
        let d1 = -1*((a*(a+b)*x) / (a * (a+1)));
        return 1/(1 + d1/(1 + recurseContinuedFractionIncompleteBeta(m+1, x, a, b, maxM)));
    }

    let d2mNum = m*(b-m)*x;
    let d2mDenom = (a+2*m-1) * (a+2*m);
    let d2m = d2mNum/d2mDenom;
    let d2m1Num = (a+m)*(a+b+m)*(x);
    let d2m1Denom = (a+2*m)*(a+2*m+1);
    let d2m1 = -1*d2m1Num/d2m1Denom;

    return (d2m/(1 + d2m1/(1 + recurseContinuedFractionIncompleteBeta(m+1, x, a, b, maxM))));
}

// https://dlmf.nist.gov/8.17#ii  8.17.22
// https://www.jstor.org/stable/2235770
function regularizedIncompleteBeta(x, a, b, maxIter=20){ // continued fraction approximation
    let firstPart = (x**a)*((1-x)**b) / (a * beta(a, b));
    let contFrac = recurseContinuedFractionIncompleteBeta(0, x, a, b, maxIter);
    return firstPart*contFrac;
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
    normalcdf,
    normalpdf,
    invNorm,
    binomcdf,
    binompdf,
    invBinom,
    tpdf,
    tcdf,
    invT,
    pearsonCorrelation,
    factorial,
    gamma,
    choose,
    erf,
    invErf,
    beta,
    regularizedIncompleteBeta,
}