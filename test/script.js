const {chiSquareGOF} = require('statlib')
// let x = [], y=[];
// for (let i = 0; i<10; i++){
//     x.push(Math.random() * 10);
//     y.push(Math.random() * 10);
// }
// console.log(geocdf(0.453, 3));

let test = [1, 2, 3, 4, 5, 6];
let test2 = [1, 3, 2, 5, 4, 6];
console.log(chiSquareGOF(test, test2));