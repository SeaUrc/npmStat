const {chiSquareTest} = require('statlib')
// let x = [], y=[];
// for (let i = 0; i<10; i++){
//     x.push(Math.random() * 10);
//     y.push(Math.random() * 10);
// }
// console.log(geocdf(0.453, 3));

let test = [[2, 2, 1],[2, 1, 2], [1, 1, 2]];
let test2 = [[4, 3, 0], [0, 3, 4], [0, 1, 0]];

console.log(chiSquareTest(test));