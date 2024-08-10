const {geocdf, geopdf} = require('statlib')
let x = [], y=[];
for (let i = 0; i<10; i++){
    x.push(Math.random() * 10);
    y.push(Math.random() * 10);
}
console.log(geocdf(0.453, 3));