const {sum, boxPlot, pearsonCorrelation, erf, normalcdf, normalpdf} = require('statlib')
let x = [], y=[];
for (let i = 0; i<10; i++){
    x.push(Math.random() * 10);
    y.push(Math.random() * 10);
}
console.log(normalpdf(-1, 0, 1));