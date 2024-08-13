const {linRegTTest} = require('statlib')

let testx = [1, 2, 3, 4, 5, 6];
let testy = [1, 3, 2, 5, 4, 6];
console.log(linRegTTest(testx, testy, 'two'));