const {ANOVA} = require('statlib')

let test1 = [1, 2, 3, 4, 5, 6];
let test2 = [1, 3, 2, 5, 4, 9];
let test3 = [3, 5, 4, 6, 8, 7]
console.log(ANOVA([test1, test2, test3]));