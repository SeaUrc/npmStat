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
    if (percentile < 0 || percentile > 1) {
        throw new Error("Percentile out of bounds");
    }

    const aSorted = a.slice().sort((a, b) => a - b);

    let N = aSorted.length;
    let rank = (N - 1) * percentile;

    const lower = Math.floor(rank);
    const rem = rank - lower;

    if (rem == 0){
        return aSorted[lower];
    }
    return aSorted[lower] + rem * (aSorted[lower+1] - aSorted[lower]);
}

function median(a) {
    if (!a.length) {
        throw new Error("Array is empty")
    }
    return quartile(a, .50);
}

function Q1(a) {
    if (!a.length) {
        throw new Error("Array is empty")
    }
    return quartile(a, .25);
}

function Q3(a) {
    if (!a.length) {
        throw new Error("Array is empty")
    }
    return quartile(a, .75);
}

function IQR(a) {
    if (!a.length) {
        throw new Error("Array is empty")
    }
    return Q3(a) - Q1(a);
}

// thresholds by the Q1 - 1.5IQR, Q3 + 1.5IQR
function hasOutliers(a){
    if (!a.length){
        throw new Error("array is empty!");
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
    return [b, outliers];
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

function skewness(a){
    if (!a.length) {
        throw new Error("Array is empty");
    }
    return 3*(mean(a) - median(a))/populationStd(a);
}

/* --- Probability Distributions --- */

/*
* Args
* uniformpdf(a, b, x)
* returns the probability density of uniform distribution at x
* */
function uniformpdf(a, b, x){
    if (a <= x && x<=b){
        return 1/(b-a);
    }
   return 0;
}

/*
* Args
* uniformpdf(a, b, x)
* returns the cumulative probability of uniform distribution at x
* */
function uniformcdf(a, b, x){
    if (x < a){
        return 0;
    }
    if (x > b){
        return 1;
    }
    return (x-a)/(b-a);
}

/*
* Args
* invUniform(a, b, p)
* returns the point x such that uniformcdf(a, b, x) = p
* */
function invUniform(a, b, p){
    if (p < 0 || p > 1){
        throw new Error("p not within interval (0, 1)");
    }
    return a + p*(b-a);
}


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
        lower = -1e99;  // Assume the lower bound is min value
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
function invNorm(p, mu = 0, sigma = 1, tail = "left") {
    // Beasley-Springer-Moro approximation
    if (p < 0 || p > 1) {
        throw new Error("Given probability must be in the range [0, 1]");
    }
    if (p === 0) {
        return -Infinity;
    }
    if (p === 1) {
        return Infinity;
    }

    let pPrime = p;
    if (tail.toLowerCase() === "right") {
        pPrime = 1 - p;
    } else if (tail.toLowerCase() === "center") {
        pPrime = 1 - (1 - p) / 2;
    } else if (tail.toLowerCase() !== "left") {
        throw new Error("Invalid tail. Available options: left, right, center.");
    }

    const a = [
        2.50662823884,
        -18.61500062529,
        41.39119773534,
        -25.44106049637
    ];
    const b = [
        -8.4735109309,
        23.08336743743,
        -21.06224101826,
        3.13082909833
    ];
    const c = [
        0.3374754822726147,
        0.9761690190917186,
        0.1607979714918209,
        0.0276438810333863,
        0.0038405729373609,
        0.0003951896511919,
        0.0000321767881768,
        0.0000002888167364,
        0.0000003960315187
    ];

    let x = pPrime - 0.5;
    let r, result;

    if (Math.abs(x) < 0.42) {
        r = x * x;
        result = x * (((a[3] * r + a[2]) * r + a[1]) * r + a[0]) /
            ((((b[3] * r + b[2]) * r + b[1]) * r + b[0]) * r + 1);
    } else {
        if (x < 0) {
            r = pPrime;
        } else {
            r = 1 - pPrime;
        }
        r = Math.log(-Math.log(r));
        result = c[0] + r * (c[1] + r * (c[2] + r * (c[3] + r * (c[4] + r * (c[5] + r * (c[6] + r * (c[7] + r * c[8])))))));
        if (x < 0) {
            result = -result;
        }
    }

    result = mu + sigma * result;

    if (tail.toLowerCase() === "center") {
        return [mu - (result - mu), result];
    }

    return result;
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
    return regularizedIncompleteBeta(1-P,N - Math.floor(K), 1 + Math.floor(K));
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
    if (N > 170){
        throw new Error("N too high");
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
        console.log(K, binompdf(N, P, K));
        cProb += binompdf(N, P, K);
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

/*
* Args
* tcdf(lower (optional), upper, degrees of freedom)
* returns the cumulative probabilty from lower to upper of the t-distribution with given degrees of freedom
* */
function tcdf(...args){
    let lower, upper, df;
    if (args.length === 2) {
        [upper, df] = args;
        lower = -1e99;  // Assume the lower bound is min value
    } else if (args.length === 3) {
        [lower, upper, df] = args;
    } else {
        throw new Error("Invalid number of arguments. Expected 2 or 3 arguments.");
    }
    let func = (t) => {return (df/(t**2 + df))};
    let cdfx = (x) => {
        if (x === 0){
            return 0.5;
        }
        if (x < 0){
            return (0.5*regularizedIncompleteBeta(func(x), df/2, 0.5));
        }
        return (1 - 0.5*regularizedIncompleteBeta(func(x), df/2, 0.5))
    };
    return cdfx(upper) - cdfx(lower);
}

/*
* Args
* invT(probability, degrees of freedom)
* returns the point x at which the cumulative probability is equal to that given from a t-distribution with given degrees of freedom
* */
function invT(p, df) {
    if (p <= 0.0 || p >= 1.0) {
        throw new Error("p must be between 0 and 1");
    }

    const q = invNorm(p);
    const q2 = q * q;
    const q3 = q * q2;
    const q4 = q2 * q2;
    const q5 = q3 * q2;
    const q6 = q3 * q3;
    const q7 = q5 * q2;
    const q8 = q4 * q4;
    const q9 = q5 * q4;

    const a = (q2 + 1) / (4 * df);
    const b = ((5 * q2 + 16) * q2 + 3) / (96 * df * df);
    const c = (((3 * q2 + 19) * q2 + 17) * q2 - 15) / (384 * df * df * df);
    const d = ((((79 * q2 + 776) * q2 + 1482) * q2 - 1920) * q2 - 945) / (92160 * df * df * df * df);
    const e = (((((27 * q2 + 339) * q2 + 930) * q2 - 1785) * q2 - 765) * q2 + 1701) / (368640 * df * df * df * df * df);
    const f = (((((((444 * q2 + 8220) * q2 + 10675) * q2 - 16699) * q2 - 32328) * q2 + 10616) * q2 + 334305) * q2 + 40455) / (10117120 * df * df * df * df * df * df);
    const g = ((((((((-7039 * q2 + 108285) * q2 + 310853) * q2 + 270650) * q2 - 468789) * q2 - 943035) * q2 + 185285) * q2 + 113165) * q2 + 91909) / (30105600 * df * df * df * df * df * df * df);
    const h = (((((((((((125385 * q2 + 2760615) * q2 + 18441175) * q2 + 15817235) * q2 - 35111603) * q2 - 45438315) * q2 + 14744335) * q2 + 49651405) * q2 - 15886035) * q2 + 20479665) * q2 + 1007767) / (132710400 * df * df * df * df * df * df * df * df));
    const i = (((((((((((-421663 * q2 + 10092330) * q2 + 62957385) * q2 + 168379960) * q2 + 20456351) * q2 - 250887169) * q2 - 210529723) * q2 + 182954229) * q2 + 106006292) * q2 + 57455870) * q2 + 2038705) / (678297600 * df * df * df * df * df * df * df * df * df));
    const j = (((((((((((-1363430 * q2 + 30574500) * q2 + 162248910) * q2 + 330813780) * q2 + 47929268) * q2 - 681253720) * q2 - 299539594) * q2 + 519202354) * q2 + 385171956) * q2 + 154064648) * q2 + 4273375) / (2033918976 * df * df * df * df * df * df * df * df * df * df));

    let t = q * (1 + a + b + c + d + e + f + g + h + i + j);

    return t;
}


/*
* Args
* chipdf(x, degrees of freedom)
* returns the probability density at x of a chi-square distribution with given degrees of freedom
* */
function chipdf(x, k){
    if (x <= 0){
        return 0;
    }
    let num = x**((k/2) - 1) * Math.exp(-x/2);
    let denom = 2**(k/2) * gamma(k/2);

    return num/denom;
}

/*
* Args
* chicdf(x, degrees of freedom)
* returns the cumulative probability up to x (left tailed) of a chi-square distribution with given degrees of freedom
* */
function chicdf(...args){
    let lower, upper, df;
    if (args.length === 2) {
        [upper, df] = args;
        lower = 0;  // Assume the lower bound is min value
    } else if (args.length === 3) {
        [lower, upper, df] = args;
        lower = Math.max(lower, 0);
    } else {
        throw new Error("Invalid number of arguments. Expected 2 or 3 arguments.");
    }

    const cdf = (x) => {
        return (lowerIncompleteGamma(df/2, x/2) / gamma(df/2));
    }
    return cdf(upper) - cdf(lower);
}

/*
* Args
* fpdf(x, df1, df2)
* returns the probability density of a F distribution with given degrees of freedom at x
* */
function fpdf(x, dfNum, dfDenom){
    if (x < 0){
        throw new Error("x cannot be <= 0!");
    }
    if (x === 0){
        return 0;
    }
    let num = Math.sqrt((((dfNum*x)**dfNum) * (dfDenom**dfDenom)) / (((dfNum**x) + dfDenom) ** (dfNum+dfDenom)));
    return num / (x * beta(dfNum*0.5, dfDenom*0.5));
}

/*
* Args
* fcdf(lower (optional), upper, df1, df2)
* returns the cumulative probability from lower to upper of a F distribution with given degrees of freedom
* */
function fcdf(...args){
    let lower, upper, df1, df2;
    if (args.length === 3) {
        [upper, df1, df2] = args;
        lower = 0;  // Assume the lower bound is min value
    } else if (args.length === 4) {
        [lower, upper, df1, df2] = args;
        lower = Math.max(lower, 0);
    } else {
        throw new Error("Invalid number of arguments. Expected 2 or 3 arguments.");
    }
    const cdf = (x) => {
        return regularizedIncompleteBeta((df1*x)/(df1*x+df2), df1/2, df2/2);
    }
    return cdf(upper)-cdf(lower);
}

/*
* Args
* poissonpdf(lambda, k)
* returns the probability mass of the poission distribution with lambda at k.
* */
function poissonpdf(lambda, k){ // should technically be called a probability mass function (PMF) but i wanted to keep the names the same
    if (lambda <=0 ){
        throw new Error("lambda must be in the interval (0, Infinity)");
    }
    if (k < 0){
        throw new Error("k cannot be < 0");
    }
    if (!Number.isInteger(k)){
        throw new Error("k must be an non-negative integer");
    }
    return (lambda**k) * Math.exp(-lambda) / factorial(k);
}

/*
* Args
* poissoncdf(lambda, k)
* returns the cumulative probability from 0 to k of the poission distribution with lambda.
* */
function poissoncdf(lambda, k){
    if (lambda <=0 ){
        throw new Error("lambda must be in the interval (0, Infinity)");
    }
    if (k < 0){
        throw new Error("k cannot be < 0");
    }
    if (!Number.isInteger(k)){
        throw new Error("k must be an non-negative integer");
    }
    return upperIncompleteGamma(Math.floor(k+1), lambda) / (factorial(Math.floor(k)));
}

/*
* Args
* geopdf(p, x)
* returns the probability that the first occurrence of success requires x independent trials, each with a probability of success of p
* */
function geopdf(p, x){
    if (p<=0 || p > 1){
        throw new Error("p must be within the range (0, 1]");
    }
    if (x < 1 || !Number.isInteger(x)){
        throw new Error("x must be a natural number");
    }
    return (1-p)**(x-1) * p;
}

/*
* Args
* geocdf(p, x)
* returns the probability that the first occurrence of success requires x or less independent trials, each with a probability of success of p
* */
function geocdf(p, x){
    if (p<=0 || p > 1){
        throw new Error("p must be within the range (0, 1]");
    }
    if (x < 1 || !Number.isInteger(x)){
        throw new Error("x must be a natural number");
    }
    return 1-((1-p)**(Math.floor(x)));
}



/* --- Correlation and Regression --- */
function pearsonCorrelation(x, y) {
    if (!x.length || !y.length) {
        throw new Error("Empty array");
    }
    if (x.length !== y.length) {
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

/*
* args
* Data: zTest(populationMean, populationStdDev, sampleList, tail)
* Stats: zTest(populationMean, populationStdDev, sampleMean, number of samples, tail)
* returns z score and p value
* */
function zTest(...args){
    if (args.length < 4 || args.length > 5){
        throw new Error("found " + args.length + " parameters but expected 4 or 5");
    }

    let stats = true;
    args.forEach((arg) => {
        if (Array.isArray(arg)){
            stats = false;
        }
    })

    let populationMean, populationStdDev, sampleMean, N, tail;
    if (stats) {
        populationMean = args[0];
        populationStdDev = args[1];
        sampleMean = args[2];
        N = args[3];
        tail = args[4];
    }else{
        populationMean = args[0];
        populationStdDev = args[1];
        sampleMean = mean(args[2]);
        N = args[2].length;
        tail = args[3];
    }

    const stdErr = populationStdDev / Math.sqrt(N);
    const zScore = (sampleMean - populationMean) / stdErr;
    let pValue;
    if (tail === 'left') {
        pValue = normalcdf(zScore, 0, 1);
    } else if (tail === 'right') {
        pValue = 1 - normalcdf(zScore, 0, 1);
    } else if (tail==='two'){
        pValue = 1 - normalcdf(-Math.abs(zScore), Math.abs(zScore), 0, 1);
    }else{
        throw new Error("unknown tail type. Use 'left', 'right', or 'two'");
    }
    return {zScore, pValue};
}

/*
* args
* Data: twoSampleZTest(sample1, sample2, tail)
* Stats: twoSampleZTest(sampleStdDev1, sampleStdDev2, mean1, sampleSize1, mean2, sampleSize2, tail)
* returns z score and p-value
* */
function twoSampleZTest(...args){
    if (args.length !== 3 && args.length !== 7){
        throw new Error("should be 3 or 7 arguments");
    }
    let mean1, mean2, variance1, variance2, sampleSize1, sampleSize2, tail;
    if (args.length === 3){
        let sample1 = args[0];
        let sample2 = args[1];
        tail = args[2];
        mean1 = mean(sample1);
        mean2 = mean(sample2);
        variance1 = sampleVariance(sample1);
        variance2 = sampleVariance(sample2);
        sampleSize1 = sample1.length;
        sampleSize2 = sample2.length;
    }else{
        mean1 = args[2];
        mean2 = args[4];
        variance1 = args[0]**2;
        variance2 = args[1]**2;
        sampleSize1 = args[3];
        sampleSize2 = args[5];
        tail = args[6];
    }

    const stderr =  Math.sqrt((variance1/sampleSize1) + (variance2/sampleSize2));;
    const zScore = (mean1-mean2) / stderr;

    let pValue;
    if (tail === 'left'){
        pValue = normalcdf(zScore, 0, 1);
    }else if (tail === 'right'){
        pValue = 1-normalcdf(zScore, 0, 1);
    }else if (tail === 'two'){
        pValue = 1-normalcdf(-Math.abs(zScore), Math.abs(zScore), 0, 1);
    }else{
        throw new Error("unknown tail type. Use 'left', 'right', or 'two'");
    }
    return {zScore, pValue};
}

/*
* args
* onePropZTest(hypothesized population proportion, success, trials, tail)
* returns z score and p-value
* */
function onePropZTest(p0, x, n, tail = "two") {
    const stdErr = Math.sqrt((p0 * (1 - p0)) / n);

    const zScore = (x/n - p0) / stdErr;

    let pValue;
    if (tail === "left") {
        pValue = normalcdf(zScore, 0, 1);
    } else if (tail === "right") {
        pValue = 1 - normalcdf(zScore, 0, 1);
    } else if (tail === "two") {
        pValue = 1 - normalcdf(-Math.abs(zScore), Math.abs(zScore), 0, 1);
    } else {
        throw new Error("Invalid tail type specified. Use 'two-tailed', 'left-tailed', or 'right-tailed'.");
    }

    return {zScore, pValue}
}

/*
* args
* twoProportionZTest(success1, sampleSize1, success2, sampleSize2, tail)
* returns z score and p-value
* */
function twoProportionZTest(x1, n1, x2, n2, tail = "two") {
    const pPooled = (x1+x2)/(n1+n2);
    const stdErr = Math.sqrt(pPooled * (1-pPooled) * (1/n1 + 1/n2));
    const zScore = (x1/n1 - x2/n2) / stdErr;

    let pValue;
    if (tail === "left") {
        pValue = normalcdf(zScore, 0, 1);
    } else if (tail === "right") {
        pValue = 1 - normalcdf(zScore, 0, 1);
    } else if (tail === "two") {
        pValue = 1 - normalcdf(-Math.abs(zScore), Math.abs(zScore), 0, 1);
    } else {
        throw new Error("Invalid tail type specified. Use 'two-tailed', 'left-tailed', or 'right-tailed'.");
    }

    return {zScore, pValue}
}

/*
* args
* Data: tTest(population mean, sample list, tail)
* Stats: tTest(population mean, sample mean, sample standard deviation, sample size, tail)
* returns z score and p-value
* */
function tTest(...args){
    if (args.length !== 3 && args.length !== 5){
        throw new Error("found " + args.length + " parameters but expected 3 or 5");
    }

    let populationMean, sampleMean, sampleStdDev, N, tail;
    if (args.length === 5) {
        populationMean = args[0];
        sampleMean = args[1];
        sampleStdDev = args[2];
        N = args[3];
        tail = args[4];
    }else{
        populationMean = args[0];
        sampleMean = mean(args[1]);
        sampleStdDev = sampleStd(args[1]);
        N = args[1].length;
        tail = args[2];
    }


    const testStat = (sampleMean - populationMean) / (sampleStdDev/(Math.sqrt(N)));
    let pValue;

    if (tail === 'left') {
        pValue = tcdf(testStat, N-1);
    } else if (tail === 'right') {
        pValue = 1-tcdf(testStat, N-1);
    } else if (tail==='two'){
        pValue = 1 - tcdf(-Math.abs(testStat), Math.abs(testStat), N-1);
    }else{
        throw new Error("unknown tail type. Use 'left', 'right', or 'two'");
    }
    return {testStat, pValue};
}

/*
* args
* Data: twoSampleTTest(sample list 1, sample list 2, pooling, tail)
* Stats: twoSampleTTest(sample mean 1, sample mean 2, sample standard deviation 1, sample standard deviation 2, sample size 1, sample size2, pooling, tail)
* returns z score and p-value
* */
function twoSampleTTest(...args){
    if (args.length !== 4 && args.length !== 8){
        throw new Error("found " + args.length + " parameters but expected 4 or 8");
    }

    let mean1, mean2, variance1, variance2, n1, n2, pooling, tail;
    if (args.length === 8) {
        mean1 = args[0];
        variance1 = args[1]**2;
        n1 = args[2];
        mean2 = args[3];
        variance2 = args[4]**2;
        n2 = args[5];
        pooling = args[6];
        tail = args[7];
    }else{
        mean1 = mean(args[0]);
        mean2 = mean(args[1]);
        variance1 = sampleVariance(args[0]);
        variance2 = sampleVariance(args[1]);
        n1 = args[0].length;
        n2 = args[1].length;
        pooling = args[2];
        tail = args[3];
    }


    let testStat;
    let df;
    if (pooling){
        const SDpooled = Math.sqrt(((n1-1) * variance1 + (n2-1) * variance2 ) / (n1+n2-2));
        const SEpooled = SDpooled * Math.sqrt(1/n1 + 1/n2);
        testStat = (mean1 - mean2) / SEpooled;
        df = n1+n2-2;
    }else{
        const SDunequal = Math.sqrt((variance1/n1) + (variance2/n2));
        testStat = (mean1 - mean2) / SDunequal;
        // Welch-Satterthwaite equation.
        const dfNum = (variance1/n1 + variance2/n2)**2;
        const dfDenom = ( (variance1**2) / ((n1**2) * (n1-1)) ) + ( (variance2**2) / ((n2**2) * (n2-1)));
        df = dfNum/dfDenom;
    }


    // const testStat = (sampleMean - populationMean) / (sampleStdDev/(Math.sqrt(N)));
    let pValue;

    if (tail === 'left') {
        pValue = tcdf(testStat, df);
    } else if (tail === 'right') {
        pValue = 1-tcdf(testStat, df);
    } else if (tail==='two'){
        pValue = 1 - tcdf(-Math.abs(testStat), Math.abs(testStat), df);
    }else{
        throw new Error("unknown tail type. Use 'left', 'right', or 'two'");
    }
    return {testStat, pValue};
}

/*
* args
* chiSquareTest(observed)
* returns expected, chi square statistics, p value
* */
function chiSquareTest(A){
    if (!A.length || !A[0].length){
        throw new Error("size of A must be greater than 1");
    }
    const numRows = A.length;
    const numCols = A[0].length;

    const rowTot = Array(numRows).fill(0);
    const colTot = Array(numCols).fill(0);
    let grandTot = 0;

    for (let i=0; i<numRows; i++){
        for (let j=0; j<numCols; j++){
            rowTot[i] += A[i][j];
            colTot[j] += A[i][j];
            grandTot += A[i][j];
        }
    }
    const B = Array.from({length: numRows}, () => Array(numCols).fill(0));
    for (let i=0; i<numRows; i++){
        for (let j=0; j<numCols; j++){
            B[i][j] = (rowTot[i] * colTot[j]) / grandTot;
        }
    }

    let chiSquare = 0;
    for (let i=0; i<numRows; i++){
        for (let j=0; j<numCols; j++){
            chiSquare += ((A[i][j] - B[i][j])**2) / B[i][j];
        }
    }

    const df = (numRows-1) * (numCols-1);

    const pValue = 1-chicdf(chiSquare, df);
    return {B, chiSquare, pValue};
}

/*
* args
* chiSquareGOF(observed, expected)
* returns chi square statistics, p value
* */
function chiSquareGOF(A, B){
    if (!A.length !== !B.length){
        throw new Error("lengths of observed and expected must be the same");
    }
    if (A.length === 0){
        throw new Error("lengths of observed and expected must be greater than 0")
    }

    const N = A.length;
    let chiSquare = 0;
    for (let i=0; i<N; i++){
        chiSquare += ((A[i] - B[i])**2) / B[i];
    }

    const df = N-1;

    const pValue = 1 - chicdf(chiSquare, df);
    return {chiSquare, pValue};
}


/*
* args
* linRegTTest(xlist, ylist, tail)
* returns slope, intercept, test statistics, p value
* */
function linRegTTest(x, y, tail='two'){
    if (x.length !== y.length){
        throw new Error("x and y are not the same length");
    }
    const N = x.length;

    const meanX = mean(x);
    const meanY = mean(y);

    let SSxy = 0;
    let SSxx = 0;
    for (let i=0; i<N; i++){
        SSxy += (x[i] - meanX) * (y[i] - meanY);
        SSxx += (x[i] - meanX)**2;
    }

    const b1 = SSxy/SSxx;
    const b0 = meanY - b1*meanX;

    const yPred = x.map((i) => {return (b0 + b1*i)});

    let resSS = 0;
    for (let i=0; i<N; i++){
        resSS += (y[i] - yPred[i]) ** 2;
    }

    const s = Math.sqrt(resSS / (N-2));
    const SEb1 = s / Math.sqrt(SSxx);

    const tStat = b1 / SEb1;
    const df = N-2;

    let pVal;
    if (tail === 'two'){
        pVal = 1 - tcdf(-Math.abs(tStat), Math.abs(tStat), df);
    }else if (tail === 'left'){
        pVal = tcdf(tStat, df);
    }else if (tail === 'right'){
        pVal = 1 - tcdf(tStat, df);
    }else{
        throw new Error("unsupported tail type. Use 'left', 'right', or 'two'");
    }

    return {
        b1,
        b0,
        tStat,
        pVal
    };
}

/*
* args
* Data: twoSampleFTest(sample1, sample2, tail);
* Stats: twoSampleFTest(standard deviation 1, sample size 1, standard deviation 2, sample size 2, tail)
* returns test statistics and p value
* */
function twoSampleFTest(...args){
    if (args.length !== 3 && args.length !== 5){
        throw new Error("found " + args.length + " parameters but expected 3 or 5");
    }

    let var1, var2, n1, n2, tail;
    if (args.length === 5){
        var1 = args[0]**2;
        n1 = args[1];
        var2 = args[2]**2;
        n2 = args[3];
        tail = args[4];
    }else{
        let sample1 = args[0];
        let sample2 = args[1];
        tail = args[2];

        var1 = sampleVariance(sample1);
        var2 = sampleVariance(sample2);
        n1 = sample1.length;
        n2 = sample2.length;
    }

    const F = var1/var2;
    const df1 = n1-1;
    const df2 = n2-1;

    let pVal;
    if (tail === 'left'){
        pVal = fcdf(F, df1, df2);
    }else if (tail === 'right'){
        pVal = 1-fcdf(F, df1, df2);
    }else if (tail === 'two'){
        if (F > 1){
            pVal = 2 * (1 - fcdf(F, df1, df2))
        }else{
            pVal = 2 * fcdf(F, df1, df2);
        }
    }else{
        throw new Error("Unsupported tail type. Use 'left', 'right', or 'two'")
    }

    return {F, pVal}
}


/*
* args
* ANOVA(samples);
* returns test statistics and p value
* */
function ANOVA(samples){
    const K = samples.length;
    let N = 0;
    for (let i=0; i<K; i++){
        N += samples[i].length;
    }

    const means = samples.map((sample) => mean(sample));
    let totMean = 0;

    samples.forEach((sample) => {
        sample.forEach((s) => {
            totMean += s;
        })
    })
    totMean /= N;

    let SSB = 0;
    let SSW = 0;
    samples.forEach((sample, i) => {
        SSB += sample.length * ((means[i] - totMean) ** 2);
        sample.forEach((s) => {
            SSW += (s - means[i]) ** 2;
        })
    })

    const dfb = K-1;
    const dfw = N-K;

    const MSB = SSB / dfb;
    const MSW = SSW / dfw;

    const F = MSB/MSW;

    const pVal = 1 - fcdf(F, dfb, dfw);

    return {F, pVal};
}

/* --- Confidence Intervals --- */

/*
* args
* Data: zInterval(population standard deviation, sample list, confidence level)
* Stats: zInterval(population standard deviation, sample mean, sample size, confidence level)
* */
function zInterval(...args){
    if (args.length !== 3 && args.length !== 4){
        throw new Error(`found ${args.length} args but only 3 or 4 are expected`)
    }
    let populationStdDev, sampleMean, N, confidence;
    if (args.length === 4){
        populationStdDev = args[0];
        sampleMean = args[1];
        N = args[2];
        confidence = args[3];
    }else{
        populationStdDev = args[0];
        sampleMean = mean(args[1]);
        N = args[1].length;
        confidence = args[2];
    }
    
    const zCrit = invNorm(confidence, 0, 1, 'center')[1]; // only care abt positive z critical
    const moe = zCrit * (populationStdDev / Math.sqrt(N));
    return [sampleMean - moe, sampleMean + moe];
}



/* --- RNG --- */
function randomUniform(min = 0, max = 1){
    return Math.random() * (max - min) + min;
}

function randomNormal(mean = 0, stdDev = 1) {
    let u1 = Math.random();
    let u2 = Math.random();
    let z = Math.sqrt(-2.0 * Math.log(u1)) * Math.cos(2.0 * Math.PI * u2);
    return z * stdDev + mean;
}

function randomBinomial(n, p) {
    if (n*p > 500 && n*(1-p) > 500){ // when n becomes large
        return randomNormal(n*p, Math.sqrt(n*p*(1-p)));
    }
    let count = 0;
    for (let i = 0; i < n; i++) {
        if (Math.random() < p) {
            count++;
        }
    }
    return count;
}

function randomPoisson(lambda) {
    let L = Math.exp(-lambda);
    let k = 0;
    let p = 1;
    do {
        k++;
        p *= Math.random();
    } while (p > L);
    return k - 1;
}

/* --- Matrix Operations --- */

function matmul(A, B){ // N1,N2 x N2xN3 = N1,N3
    let result = Array(A.length).fill().map(() => Array(B[0].length).fill(0)); // N1,N3
    for (let i = 0; i < A.length; i++) {
        for (let j = 0; j < B[0].length; j++) {
            for (let k = 0; k < B.length; k++) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return result;
}

function transpose(A){
    return A[0].map((_, colIndex) => A.map(row => row[colIndex]));
}

function inverseMatrix(A) { // matrix x answer = identity
    const size = A.length;
    const identity = Array(size).fill().map((_, i) => Array(size).fill(0).map((_, j) => (i === j ? 1 : 0)));

    for (let i = 0; i < size; i++) {
        let diag = A[i][i];
        for (let j = 0; j < size; j++) {
            A[i][j] /= diag;
            identity[i][j] /= diag;
        }
        for (let k = 0; k < size; k++) {
            if (k !== i) {
                let factor = A[k][i];
                for (let j = 0; j < size; j++) {
                    A[k][j] -= factor * A[i][j];
                    identity[k][j] -= factor * identity[i][j];
                }
            }
        }
    }
    return identity;
}

/* --- Misc --- */
function factorial(x) {
    if (!Number.isInteger(x)) {
        throw new Error("Non integer factorial. If this was intentional, use the gamma function")
    }
    if (x === 0){
        return 1;
    }
    if (x === 1) {
        return 1;
    }
    return x * factorial(x - 1);
}

function choose(N, K){
    if (!Number.isInteger(N) || !Number.isInteger(K)){
        throw new Error("Non integer in combinatoric");
    }
    return ((factorial(N)) / (factorial(K) * factorial(N-K)));
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

function lowerIncompleteGamma(s, z, kLimit=50){
    let first = z**s * Math.exp(-z);
    let powerSeries = 0;
    for (let k = 0; k<=kLimit; k++){
        powerSeries += (z**k) / pochhammer(s, k+1);
    }
    return first * powerSeries;
}

// lower + upper incomplete gamma = gamma
function upperIncompleteGamma(s, z, kLimit=50){
    return gamma(s) - lowerIncompleteGamma(s, z, kLimit);
}

// Abramowitz and Stegun approximation
function erf(x) {
    const a1 =  0.254829592;
    const a2 = -0.284496736;
    const a3 =  1.421413741;
    const a4 = -1.453152027;
    const a5 =  1.061405429;
    const p =  0.3275911;

    const sign = x < 0 ? -1 : 1;
    x = Math.abs(x);

    const t = 1.0 / (1.0 + p * x);
    const y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.exp(-x * x);

    return sign * y;
}

// using a taylor series approximation
function invErf(s) {
    // copied from wolfram alpha ~ 60 terms.
    return ((Math.sqrt(Math.PI)*s)/2 + 1/24*Math.PI**(3/2)*s**3 + 7/960*Math.PI**(5/2)*s**5 + (127*Math.PI**(7/2)*s**7)/80640 + (4369*Math.PI**(9/2)*s**9)/11612160 + (34807*Math.PI**(11/2)*s**11)/364953600 + (20036983*Math.PI**(13/2)*s**13)/797058662400 + (2280356863*Math.PI**(15/2)*s**15)/334764638208000 + (49020204823*Math.PI**(17/2)*s**17)/26015994740736000 + (65967241200001*Math.PI**(19/2)*s**19)/124564582818643968000 + (15773461423793767*Math.PI**(21/2)*s**21)/104634249567660933120000 + (655889589032992201*Math.PI**(23/2)*s**23)/15127122937496123473920000 + (94020690191035873697*Math.PI**(25/2)*s**25)/7474578392645143363584000000 + (655782249799531714375489*Math.PI**(27/2)*s**27)/178403237075654281802022912000000 + (44737200694996264619809969*Math.PI**(29/2)*s**29)/41389551001551793378069315584000000 + (10129509912509255673830968079*Math.PI**(31/2)*s**31)/31699526708247314693086028759040000000 + (108026349476762041127839800617281*Math.PI**(33/2)*s**33)/1138139806932911586740560776564572160000000 + (10954814567103825758202995557819063*Math.PI**(35/2)*s**35)/386967534357189939491790664031954534400000000 + (61154674195324330125295778531172438727*Math.PI**(37/2)*s**37)/7216170580692877991642912302867888157491200000000 + (54441029530574028687402753586278549396607*Math.PI**(39/2)*s**39)/21388729601173690367229592065700420498803916800000000 + (452015832786609665624579410056180824562551*Math.PI**(41/2)*s**41)/589538093208821045416076151054599825513250816000000000 + (2551405765475004343830620568825540664310892263*Math.PI**(43/2)*s**43)/11017390414250484013439181732847675382639546597376000000000 + (70358041406630998834159902148730577164631303295543*Math.PI**(45/2)*s**45)/1003463918929934083944040672227766273850809904089006080000000000 + (775752883029173334450858052496704319194646607263417*Math.PI**(47/2)*s**47)/36461999877756596461966654342124885446478168279671111680000000000 + (132034545522738294934559794712527229683368402215775110881*Math.PI**(49/2)*s**49)/20410552443571076541093845901017299875687763126921176211783680000000000 + (205123592348209083599687465864398200639531076044966484944001*Math.PI**(51/2)*s**51)/104093817462212490359578614095188229366007591947297998680096768000000000000 + (49286712796548457688918350570407306376864714936222866274867681*Math.PI**(53/2)*s**53)/81966445978816463837428188698953931466490549544786652674956197888000000000000 + (626335302110212978311199487037804496444985797387366410523685237887*Math.PI**(55/2)*s**55)/3408164823799188566360264086102504470376677050072229018224678708183040000000000000 + (71991822498077532856239242609441182504985736043177831717229132014817*Math.PI**(57/2)*s**57)/1279866131478471753390819172098728737581453311038888826608608757237678080000000000000 + (11830685611118801673299916015156726691775984754422156198050789985296297*Math.PI**(59/2)*s**59)/686220574494272956513894143020013323318557280335265609208564109157093867520000000000000);
}

function beta(z1, z2){
    return (gamma(z1)*gamma(z2))/gamma(z1+z2);
}

function recurseContinuedFractionIncompleteBeta(m, x, a, b, maxM){
    if (m >= maxM){
        return 1;
    }
    if (m === 0){
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
function regularizedIncompleteBeta(x, a, b, maxIter=50){ // continued fraction approximation
    let firstPart = (x**a)*((1-x)**b) / (a * beta(a, b));
    let contFrac = recurseContinuedFractionIncompleteBeta(0, x, a, b, maxIter);
    return firstPart*contFrac;
}

function pochhammer(x, n){
    let res = 1;
    for (let i = x; i<=(x+n-1); i++){
        res *= i;
    }
    return res;
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
    hasOutliers,
    sum,
    mean,
    sampleVariance,
    populationVariance,
    sampleStd,
    populationStd,
    coefficientOfVariation,
    skewness,
    uniformpdf,
    uniformcdf,
    invUniform,
    normalpdf,
    normalcdf,
    invNorm,
    binompdf,
    binomcdf,
    invBinom,
    tpdf,
    tcdf,
    invT,
    chipdf,
    chicdf,
    fpdf,
    fcdf,
    poissonpdf,
    poissoncdf,
    geopdf,
    geocdf,
    pearsonCorrelation,
    zTest,
    twoSampleZTest,
    onePropZTest,
    twoProportionZTest,
    tTest,
    twoSampleTTest,
    chiSquareTest,
    chiSquareGOF,
    linRegTTest,
    twoSampleFTest,
    ANOVA,
    zInterval,
    randomUniform,
    randomNormal,
    randomBinomial,
    randomPoisson,
    matmul,
    transpose,
    inverseMatrix,
    factorial,
    choose,
    gamma,
    lowerIncompleteGamma,
    upperIncompleteGamma,
    erf,
    invErf,
    beta,
    regularizedIncompleteBeta,
    pochhammer,
}