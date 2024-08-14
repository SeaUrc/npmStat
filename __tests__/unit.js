const {min, max, range, quartile, median, Q1, Q3, IQR, hasOutliers, sum, sampleVariance} = require("../index");

const accuracy = 5;

describe('Descriptive funcs', () => {
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
        expect(range([1.4, 3.23, 5.142, 6.1253, 2.06912, 3.9623, 4.2535])).toBeCloseTo(4.7253, accuracy);
    })

    test('quartile()', () => {
        expect(quartile([1, 2, 3, 4, 5], .5)).toBe(3);
        expect(quartile([4, 2, 1, 8, 16], .3)).toBe(2.4);
    })

    test('median()', () => {
        expect(median([5, 1, 3, 2, 4])).toBe(3);
        expect(median([1, 5, 23.5, 3, 100, 50])).toBeCloseTo(14.25, accuracy);
    })

    test('Q1()', () => {
        expect(Q1([5, 1, 3, 2, 4])).toBe(2);
        expect(Q1([10, 5, 3, 8, 5, 4])).toBeCloseTo(4.25, accuracy);
    })

    test('Q3()', () => {
        expect(Q3([5, 1, 3, 2, 4])).toBe(4);
        expect(Q3([10, 5, 3, 8, 5, 4])).toBeCloseTo(7.25, accuracy);
    })

    test('IQR()', () => {
        expect(IQR([1, 6, 3, 4, 56, 2])).toBeCloseTo(3.25, accuracy);
    })

    test('hasOutliers()', () => {
        expect(hasOutliers([1, 2, 3, 4, 5, 6, 7, 8, 9])).toStrictEqual([[1, 2, 3, 4, 5, 6, 7, 8, 9], []]);
        expect(hasOutliers([4, 5, 1, 3, 2, 3,5, 10, 100])).toStrictEqual([[4, 5, 1, 3, 2, 3,5], [10, 100]])
    })

    test('sum', () => {
        expect(sum([1, 2, 3, 4, 5, 6])).toBeCloseTo(21, accuracy);
        expect(sum([-103, 24.43, 651.1235, 16.342, -523.2353])).toBeCloseTo(65.6602, accuracy);
    })

    test('sample variance', () => {
        expect(sampleVariance([1, 5, 3, 4, 6, 3])).toBeCloseTo(3.06666642, accuracy);
        expect(sampleVariance([1.432, -1.532, 4.423, 1.534, 3.2423, -2.124, -1.324, 0.423])).toBeCloseTo(5.52248332, accuracy)
    })

    // test('sample standard deviation', () => {})
})