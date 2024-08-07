function min(a, comp = (u, v) => {return u<v}){
    if (!a.length){
        throw new Error("Array is empty");
    }
    let mn = a[0];
    for (let i=0; i<a.length; i++){
        if (comp(a[i], mn)){
            mn = a[i];
        }
    }
    return mn
}

function max(a, comp = (u, v) => {return u>v}){
    if (!a.length){
        throw new Error("Array is empty");
    }
    let mn = a[0];
    for (let i=0; i<a.length; i++){
        if (comp(a[i], mn)){
            mn = a[i];
        }
    }
    return mn
}

function range(a){
    if (!a.length){
        throw new Error("Array is empty");
    }
    return max(a) - min(a);
}

function quartile(a, percentile) {
    if (!a.length) {
        throw new Error("Array is empty");
    }
    if (percentile < 0 || percentile > 100) {
        throw new Error("Percentile out of bounds");
    }

    a.sort((a, b) => a - b);

    let N = a.length;
    let rank = (percentile / 100) * (N); // 0 -> N

    if (rank < 1){
        return a[0];
    }
    rank--;

    let lowerIndex = Math.floor(rank);
    let upperIndex = Math.ceil(rank);
    let lowerValue = a[lowerIndex];
    let upperValue = a[upperIndex];
    let frac = rank - Math.floor(rank);
    return lowerValue + (upperValue - lowerValue) * frac;
}

function median(a){
    if (!a.length){
        throw new Error("Array is empty")
    }
    return quartile(a, 50);
}

function Q1(a){
    if (!a.length){
        throw new Error("Array is empty")
    }
    return quartile(a, 25);
}

function Q3(a){
    if (!a.length){
        throw new Error("Array is empty")
    }
    return quartile(a, 75);
}

function IQR(a){
    if (!a.length){
        throw new Error("Array is empty")
    }
    return Q3(a)-Q1(a);
}


module.exports = {
    min,
    max,
    range,
    quartile,
    median,
    Q1,
    Q3,
    IQR
}