const {sum, boxPlot} = require('statlib')
let x = [];
for (let i = 0; i<10; i++){
    x.push(Math.random() * 10);
}
let [x1, x2, x3, x4, x5, x6] = boxPlot(x, true);
console.log(x3, x6)