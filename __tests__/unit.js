const {min, max, range, quartile} = require("../index");

describe('Descriptive funcs', () => {
    test('min() should correctly give the min of array', () => {
        expect(min([4, 3.3, 1, 2.4, 3.9, 4])).toBe(1);
        expect(min([10.4, 3, 5.3, -3.34, 4, 1.2, 6], (u, v) => {return u > v})).toBe(10.4);
    })

    test('max() should correctly give the max of an array', () => {
        expect(max([7.34, 3.12, 4.95, 2.13, 2.01, 4.10])).toBe(7.34);
        expect(max([1.34, -3.12, 2.95, 3.13, -6.01, 4.10], (u, v) => u<v)).toBe(-6.01);
    })

    test('range should give range of an array', () => {
        expect(range([1, 3, 5, 6, 2, 3, 4])).toBe(5);
        expect(range([1.4, 3.23, 5.142, 6.1253, 2.06912, 3.9623, 4.2535])).toBeCloseTo(4.7253, 5);
    })

    test('quartile should give the given quartile', () => {
        expect(quartile([1, 2, 3, 4, 5], 20)).toBe(1);
        // expect(quartile([1, 2, 4, 8, 16], 50)).toBe(4);
    })
})