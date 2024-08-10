const {geocdf, geopdf, randomUniform, randomNormal, randomBinomial} = require('statlib')
const {binomcdf} = require("../index");
// let x = [], y=[];
// for (let i = 0; i<10; i++){
//     x.push(Math.random() * 10);
//     y.push(Math.random() * 10);
// }
// console.log(geocdf(0.453, 3));

let test = []
let mx = 10000000;
let buckets = 10;
for (let i =0; i<buckets; i++){
    test.push(0);
}
for (let i=0; i<mx; i++){

}

console.log(test);