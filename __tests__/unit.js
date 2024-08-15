const {min, max, range, quartile, median, Q1, Q3, IQR, hasOutliers, sum, sampleVariance, sampleStd, populationVariance,
    populationStd, coefficientOfVariation, skewness, uniform, uniformpdf, uniformcdf, invUniform, normalpdf, normalcdf,
    invNorm, invErf, binompdf, binomcdf, invBinom, tpdf, tcdf, invT
} = require("../index");

const accuracyDescriptive = 7;

describe('Descriptive stats', () => {
    test('min()', () => {
        expect(min([4, 3.3, 1, 2.4, 3.9, 4])).toBe(1);
        expect(min([10.4, 3, 5.3, -3.34, 4, 1.2, 6], (u, v) => {return u > v})).toBe(10.4);
    })

    test('max()', () => {
        expect(max([7.34, 3.12, 4.95, 2.13, 2.01, 4.10])).toBe(7.34);
        expect(max([1.34, -3.12, 2.95, 3.13, -6.01, 4.10], (u, v) => u<v)).toBe(-6.01);
    })

    test('range()', () => {
        expect(range([1, 3, 5, 6, 2, 3, 4])).toBe(5);
        expect(range([1.4, 3.23, 5.142, 6.1253, 2.06912, 3.9623, 4.2535])).toBeCloseTo(4.7253, accuracyDescriptive);
    })

    test('quartile()', () => {
        expect(quartile([1, 2, 3, 4, 5], .5)).toBe(3);
        expect(quartile([4, 2, 1, 8, 16], .3)).toBe(2.4);
    })

    test('median()', () => {
        expect(median([5, 1, 3, 2, 4])).toBe(3);
        expect(median([1, 5, 23.5, 3, 100, 50])).toBeCloseTo(14.25, accuracyDescriptive);
    })

    test('Q1()', () => {
        expect(Q1([5, 1, 3, 2, 4])).toBe(2);
        expect(Q1([10, 5, 3, 8, 5, 4])).toBeCloseTo(4.25, accuracyDescriptive);
    })

    test('Q3()', () => {
        expect(Q3([5, 1, 3, 2, 4])).toBe(4);
        expect(Q3([10, 5, 3, 8, 5, 4])).toBeCloseTo(7.25, accuracyDescriptive);
    })

    test('IQR()', () => {
        expect(IQR([1, 6, 3, 4, 56, 2])).toBeCloseTo(3.25, accuracyDescriptive);
    })

    test('hasOutliers()', () => {
        expect(hasOutliers([1, 2, 3, 4, 5, 6, 7, 8, 9])).toStrictEqual([[1, 2, 3, 4, 5, 6, 7, 8, 9], []]);
        expect(hasOutliers([4, 5, 1, 3, 2, 3,5, 10, 100])).toStrictEqual([[4, 5, 1, 3, 2, 3,5], [10, 100]]);
    })

    test('sum', () => {
        expect(sum([1, 2, 3, 4, 5, 6])).toBeCloseTo(21, accuracyDescriptive);
        expect(sum([-103, 24.43, 651.1235, 16.342, -523.2353])).toBeCloseTo(65.6602, accuracyDescriptive);
    })

    test('sample variance', () => {
        expect(sampleVariance([1, 5, 3, 4, 6, 3])).toBeCloseTo(3.06666666666, accuracyDescriptive);
        expect(sampleVariance([1.432, -1.532, 4.423, 1.534, 3.2423, -2.124, -1.324, 0.423])).toBeCloseTo(5.52248332, accuracyDescriptive);
    })

    test('population variance', () => {
        expect(populationVariance([1.432, -1.532, 4.423, 1.534, 3.2423, -2.124, -1.324, 0.423])).toBeCloseTo(4.832172905, accuracyDescriptive);
        expect(populationVariance([1, 5, 6, 3, 7, 19])).toBeCloseTo(33.472222222, accuracyDescriptive);
    })

    test('population standard deviation', () => {
        expect(populationStd([1.432, -1.532, 4.423, 1.534, 3.2423, -2.124, -1.324, 0.423])).toBeCloseTo(2.198220395, accuracyDescriptive);
        expect(populationStd([1, 5, 6, 3, 7, 19])).toBeCloseTo(5.785518319, accuracyDescriptive);
    })

    test('sample standard deviation', () => {
        expect(sampleStd([1.432, -1.532, 4.423, 1.534, 3.2423, -2.124, -1.324, 0.423])).toBeCloseTo(2.349996451, accuracyDescriptive);
        expect(sampleStd([1, 5, 3, 4, 6, 3])).toBeCloseTo(1.75119007154, accuracyDescriptive);
    })

    test('coefficient of variation', () => {
        expect(coefficientOfVariation([1, 4, 3, 2, 1.2, 4.3, -1.3, 4, 2, 3, .24, 5.3, 4.2])).toBeCloseTo(0.717286725, accuracyDescriptive);
    })

    test('skewness', () => {
        expect(skewness([1, 2, 4.3, 6.2, 1.4, 3.4, 4.2, 1.342, 5.243, 5.1251, 1.1241, 5.123])).toBeCloseTo(-0.7063589127401594, accuracyDescriptive);
    })
})


const accuracyProb = 7;

describe('Probability Distributions', () => {
    test('uniform pdf', () => {
        expect(uniformpdf(1, 2, .1)).toBe(0);
        expect(uniformpdf(4, 6, 5)).toBe(.5);
        expect(uniformpdf(2, 3, 5)).toBe(0);
    })

    test('uniform cdf', () => {
        expect(uniformcdf(1, 2, .1)).toBe(0);
        expect(uniformcdf(4, 6, 5.43)).toBeCloseTo(0.715, accuracyProb);
        expect(uniformcdf(2, 3, 5)).toBe(1);
    })

    test("inverse uniform", () => {
        expect(invUniform( 3, 4, 0)).toBe(3);
        expect(invUniform(6, 9, .413)).toBeCloseTo(7.239, accuracyProb);
        expect(invUniform(2, 3, 1)).toBe(3);
    })

    test("normal pdf", () => { // 7 decimal
        expect(normalpdf( 0, 0, 1)).toBeCloseTo(0.3989422804, accuracyProb);
        expect(normalpdf(3.23, 6, 9.26)).toBeCloseTo(0.0411972499, accuracyProb);
    })

    test("normal cdf", () => {
        expect(normalcdf( 1, 0, 1)).toBeCloseTo(0.8413447404, accuracyProb);
        expect(normalcdf(1, 3, 0, 1)).toBeCloseTo(0.1573052923, accuracyProb);
        expect(normalcdf(-12.3, 313.5, 6.54, 100.32)).toBeCloseTo(0.5733754335, accuracyProb);
    })

    test("inverse normal", () => {
        expect(invNorm( .95, 0, 1, 'left')).toBeCloseTo(1.644853626, accuracyProb);
        expect(invNorm(.95, 1.43, 3.42, 'right')).toBeCloseTo(-4.195399401, accuracyProb);

        const [l, r] = invNorm(.532, 43.24, 6.54, 'center');
        expect(l).toBeCloseTo(38.49367984, accuracyProb);
        expect(r).toBeCloseTo(47.98632016, accuracyProb);
    })

    test("binomial pdf", () => {
        expect(binompdf(3, .4, 2)).toBeCloseTo(.288, accuracyProb);
        expect(binompdf(100, .423, 42)).toBeCloseTo(0.0804723057, accuracyProb);
    })

    test("binomial cdf", () => {
        expect(binomcdf(3, .4, 2)).toBeCloseTo(.936, accuracyProb);
        expect(binomcdf(100, .423, 39)).toBeCloseTo(0.2867591526, accuracyProb);
    })

    test('invBinom', () => {
        expect(invBinom(.5, 20, .6)).toBe(12);
        expect(invBinom(.6, 170, .423)).toBe(74);
    })

    test('tpdf', () => {
        expect(tpdf(3, 3)).toBeCloseTo(.0229720373, accuracyProb);
        expect(tpdf(2.123, 20.523)).toBeCloseTo(0.0465308675, accuracyProb);
    })

    test('tcdf', () => {
        expect(tcdf(3, 3)).toBeCloseTo(0.9711655572, accuracyProb);
        expect(tcdf(-1.23, 1.423, 4.23)).toBeCloseTo(.7466037916, accuracyProb);
        expect(tcdf(-1, -.5, 23)).toBeCloseTo(0.1470547344, accuracyProb);
    })

    test('invT', () => { // accurate to 5
        expect(invT(.513, 3)).toBeCloseTo(0.0353789245, 5);
        expect(invT(.234, 63.423)).toBeCloseTo(-0.7301285197, 5);
    })
})

