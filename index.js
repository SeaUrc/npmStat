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
function invNorm(p, mu=0, sigma=1, tail="lower") {
    if (p < 0 || p > 1) {
        throw new Error("Given probability must be in the range [0, 1]");
    }
    if (p == 0){
        return -1e99;
    }
    if (p == 1){
        return 1e99;
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
    const a = (q2 + 1) / (4 * df);
    const b = ((5 * q2 + 16) * q2 + 3) / (96 * df * df);
    const c = (((3 * q2 + 19) * q2 + 17) * q2 - 15) / (384 * df * df * df);
    const d = ((((79 * q2 + 776) * q2 + 1482) * q2 - 1920) * q2 - 945) / (92160 * df * df * df * df);

    const quantile = q * (1 + a + b + c + d);
    return quantile;
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
    if (x == 0){
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
    if (args.length != 3 && args.length != 7){
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
    if (args.length != 3 && args.length != 5){
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
    if (args.length != 4 && args.length != 8){
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
    if (!A.length != !B.length){
        throw new Error("lengths of observed and expected must be the same");
    }
    if (A.length == 0){
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
    if (x.length != y.length){
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
    if (tail == 'two'){
        pVal = 1 - tcdf(-Math.abs(tStat), Math.abs(tStat), df);
    }else if (tail == 'left'){
        pVal = tcdf(tStat, df);
    }else if (tail == 'right'){
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

// TODO 2-sample F test, ANOVA

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
    if (x == 1) {
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
    boxPlot,
    sum,
    mean,
    sampleVariance,
    populationVariance,
    sampleStd,
    populationStd,
    coefficientOfVariation,
    skewness,
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