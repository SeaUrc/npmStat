const {invErf, invNorm} = require('statlib')
let x = [], y=[];
for (let i = 0; i<10; i++){
    x.push(Math.random() * 10);
    y.push(Math.random() * 10);
}
console.log(invNorm(.9587, 6.5, 9.27, "center"));