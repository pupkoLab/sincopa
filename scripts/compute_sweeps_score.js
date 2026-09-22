/*
 * Browser implementation of compute_sweeps_score.py
 *
 * No plotting is performed here.
 *
 * Main function:
 *
 *   computeSweepsScore(msaText, homoplasyText, windowSize)
 *
 * Returns:
 *
 *   {
 *       scores,
 *       scoresText,
 *       summaryText
 *   }
 */


/*
 * Parse FASTA into:
 *
 * {
 *     headers: [...],
 *     sequences: [...]
 * }
 */
function parseSweepsMsa(msaText) {

    const headers = [];
    const sequences = [];

    let currentHeader = null;
    let currentSequence = "";

    const lines =
        msaText.split(/\r?\n/);

    for (const rawLine of lines) {

        const line =
            rawLine.trim();

        if (!line) {
            continue;
        }

        if (line.startsWith(">")) {

            if (currentHeader !== null) {
                headers.push(currentHeader);
                sequences.push(
                    currentSequence.toUpperCase()
                );
            }

            currentHeader =
                line.substring(1)
                    .trim()
                    .split(/\s+/)[0];

            currentSequence = "";
        }
        else {

            if (currentHeader === null) {
                throw new Error(
                    "Invalid FASTA: sequence before first header."
                );
            }

            currentSequence += line;
        }
    }


    if (currentHeader !== null) {

        headers.push(currentHeader);

        sequences.push(
            currentSequence.toUpperCase()
        );
    }


    if (sequences.length === 0) {
        throw new Error(
            "MSA contains no sequences."
        );
    }


    const msaLength =
        sequences[0].length;


    for (let i = 0; i < sequences.length; i++) {

        if (sequences[i].length !== msaLength) {

            throw new Error(
                "MSA sequences have different lengths."
            );
        }
    }


    return {
        headers,
        sequences,
        msaLength
    };
}


/*
 * Equivalent to load_homoplasy().
 */
function parseHomoplasy(homoplasyText) {

    const values =
        homoplasyText
            .trim()
            .split(/\s+/)
            .filter(x => x.length > 0)
            .map(Number);


    for (const value of values) {

        if (!Number.isInteger(value)) {
            throw new Error(
                "Invalid value in homoplasy file."
            );
        }
    }


    return values;
}


/*
 * Equivalent to:
 *
 * get_alignment_columns_from_dict_msa()
 */
function getAlignmentColumns(sequences, msaLength) {

    const columns = [];

    for (let j = 0; j < msaLength; j++) {

        let column = "";

        for (const sequence of sequences) {
            column += sequence[j];
        }

        columns.push(column);
    }

    return columns;
}


/*
 * Equivalent to Python's Counter.
 */
function countCharacters(text) {

    const counts = new Map();

    for (const char of text) {

        counts.set(
            char,
            (counts.get(char) || 0) + 1
        );
    }

    return counts;
}


/*
 * Equivalent to numpy.argmax() for a small JS array.
 *
 * numpy.argmax returns the FIRST maximum.
 */
function argmax(array) {

    let maxIndex = 0;
    let maxValue = array[0];

    for (let i = 1; i < array.length; i++) {

        if (array[i] > maxValue) {

            maxValue = array[i];
            maxIndex = i;
        }
    }

    return maxIndex;
}


/*
 * Equivalent to calculate_greedy_association().
 */
function calculateGreedyAssociation(col1, col2) {

    const charToNum = {
        "A": 0,
        "C": 1,
        "G": 2,
        "T": 3,
        "-": 4
    };


    const pairsOccurrences =
        Array.from(
            { length: 5 },
            () => Array(5).fill(0)
        );


    for (let i = 0; i < col1.length; i++) {

        const row =
            charToNum[col1[i]];

        const column =
            charToNum[col2[i]];


        if (row === undefined ||
            column === undefined) {

            throw new Error(
                "Unexpected character in fixed MSA. " +
                "Expected only A, C, G, T or -."
            );
        }


        pairsOccurrences[row][column]++;
    }


    const counts =
        countCharacters(col1);


    /*
     * Python:
     *
     * sorted(counts.items(),
     *        key=lambda x: x[1],
     *        reverse=True)
     *
     * JS stable sort preserves original insertion order
     * for equal counts.
     */
    const sortedCounts =
        Array.from(counts.entries())
            .sort(
                (a, b) => b[1] - a[1]
            )
            .map(
                item => item[0]
            );


    let result = 0;


    for (const char of sortedCounts) {

        const row =
            charToNum[char];

        const maxOccurrenceIndex =
            argmax(
                pairsOccurrences[row]
            );


        result +=
            pairsOccurrences[row]
                [maxOccurrenceIndex];


        /*
         * Prevent another character from being
         * assigned to the same partner character.
         */
        for (let i = 0; i < 5; i++) {

            pairsOccurrences[i]
                [maxOccurrenceIndex] = -1;
        }
    }


    return result / col1.length;
}


/*
 * Equivalent to:
 *
 * calculate_symmetric_greedy_association()
 */
function calculateSymmetricGreedyAssociation(
    col1,
    col2
) {

    return (
        calculateGreedyAssociation(
            col1,
            col2
        )
        +
        calculateGreedyAssociation(
            col2,
            col1
        )
    ) / 2;
}


/*
 * Equivalent to:
 *
 * calculate_total_edge_contribution()
 */
function calculateTotalEdgeContribution(
    edgeColumn,
    windowColumnsWithoutEdge,
    windowHomoplasyWithoutEdge
) {

    let edgeContributions = 0;


    for (
        let i = 0;
        i < windowColumnsWithoutEdge.length;
        i++
    ) {

        if (windowHomoplasyWithoutEdge[i]) {

            edgeContributions +=
                calculateSymmetricGreedyAssociation(
                    edgeColumn,
                    windowColumnsWithoutEdge[i]
                );
        }
    }


    return edgeContributions;
}


/*
 * Equivalent to get_first_window_score().
 */
function getFirstWindowScore(
    windowColumns,
    windowHomoplasy
) {

    let windowScore = 0;


    const homoplasiousColumns = [];


    for (
        let i = 0;
        i < windowColumns.length;
        i++
    ) {

        if (windowHomoplasy[i]) {

            homoplasiousColumns.push(
                windowColumns[i]
            );
        }
    }


    for (
        let i = 0;
        i < homoplasiousColumns.length - 1;
        i++
    ) {

        for (
            let j = i + 1;
            j < homoplasiousColumns.length;
            j++
        ) {

            windowScore +=
                calculateSymmetricGreedyAssociation(
                    homoplasiousColumns[i],
                    homoplasiousColumns[j]
                );
        }
    }


    return windowScore;
}


/*
 * Equivalent to get_next_window_score().
 */
function getNextWindowScore(
    previousWindowScore,
    previousWindowColumns,
    previousWindowHomoplasy,
    windowColumns,
    windowHomoplasy
) {

    let leftmostEdgeContribution = 0;
    let rightmostEdgeContribution = 0;


    if (previousWindowHomoplasy[0]) {

        const edgeColumn =
            previousWindowColumns[0];


        leftmostEdgeContribution =
            calculateTotalEdgeContribution(
                edgeColumn,
                previousWindowColumns.slice(1),
                previousWindowHomoplasy.slice(1)
            );
    }


    if (
        windowHomoplasy[
            windowHomoplasy.length - 1
        ]
    ) {

        const edgeColumn =
            windowColumns[
                windowColumns.length - 1
            ];


        rightmostEdgeContribution =
            calculateTotalEdgeContribution(
                edgeColumn,
                windowColumns.slice(0, -1),
                windowHomoplasy.slice(0, -1)
            );
    }


    return (
        previousWindowScore
        -
        leftmostEdgeContribution
        +
        rightmostEdgeContribution
    );
}


/*
 * Equivalent to get_contig_scores()
 * using symmetric_greedy.
 */
function getContigScores(
    msaColumns,
    msaLength,
    homoplasy,
    windowSize
) {

    const numberOfWindows =
        msaLength - windowSize + 1;


    if (numberOfWindows <= 0) {

        throw new Error(
            "Window size is larger than MSA length."
        );
    }


    const windowScores =
        new Array(numberOfWindows).fill(0);


    /*
     * First window.
     */
    let windowColumns =
        msaColumns.slice(
            0,
            windowSize
        );

    let windowHomoplasy =
        homoplasy.slice(
            0,
            windowSize
        );


    windowScores[0] =
        getFirstWindowScore(
            windowColumns,
            windowHomoplasy
        );


    /*
     * Remaining sliding windows.
     */
    for (
        let i = 1;
        i < numberOfWindows;
        i++
    ) {

        const previousWindowColumns =
            windowColumns;

        const previousWindowHomoplasy =
            windowHomoplasy;


        windowColumns =
            msaColumns.slice(
                i,
                i + windowSize
            );

        windowHomoplasy =
            homoplasy.slice(
                i,
                i + windowSize
            );


        windowScores[i] =
            getNextWindowScore(
                windowScores[i - 1],
                previousWindowColumns,
                previousWindowHomoplasy,
                windowColumns,
                windowHomoplasy
            );
    }


    /*
     * Python:
     *
     * normalization_factor =
     *     2 / (window_size *
     *          (window_size - 1))
     */
    const normalizationFactor =
        2 /
        (
            windowSize *
            (windowSize - 1)
        );


    return windowScores.map(
        score =>
            score * normalizationFactor
    );
}


/*
 * Hamming distance between two strings.
 *
 * Equivalent to get_pairwise_distance().
 */
function getPairwiseDistance(seq1, seq2) {

    let distance = 0;

    const length =
        Math.min(
            seq1.length,
            seq2.length
        );


    for (let i = 0; i < length; i++) {

        if (seq1[i] !== seq2[i]) {
            distance++;
        }
    }


    return distance;
}


/*
 * This reproduces the intended APD/pi calculation
 * using sequence values and their frequencies.
 */
function getAveragePairwiseDistanceAndPi(
    headers,
    sequences,
    msaLength
) {

    const numberOfSpecies = headers.length;

    const numOfPairs =
        numberOfSpecies * (numberOfSpecies - 1) / 2;

    let totalRelativePairwiseDistance = 0;
    let pi = 0;

    /*
     * Reproduce the Python code exactly:
     *
     * sequence2count[header] =
     *     sequence2count.get(sequence, 0) + 1
     *
     * Since the dictionary is keyed by headers, looking up
     * a sequence normally returns 0, so each header gets 1.
     */
    const sequence2count = new Map();

    for (let i = 0; i < headers.length; i++) {

        const header = headers[i];
        const sequence = sequences[i];

        sequence2count.set(
            header,
            (sequence2count.get(sequence) || 0) + 1
        );
    }

    const sequences2frequency = new Map();

    for (const [key, count] of sequence2count) {

        sequences2frequency.set(
            key,
            count / numberOfSpecies
        );
    }

    const keys =
        Array.from(sequence2count.keys());

    for (let i = 0; i < keys.length - 1; i++) {

        for (let j = i + 1; j < keys.length; j++) {

            const seq1 = keys[i];
            const seq2 = keys[j];

            const apd =
                getPairwiseDistance(
                    seq1,
                    seq2
                ) / msaLength;

            totalRelativePairwiseDistance += apd;

            pi +=
                2 *
                sequences2frequency.get(seq1) *
                sequences2frequency.get(seq2) *
                apd;
        }
    }

    return {
        apd:
            totalRelativePairwiseDistance /
            numOfPairs,

        pi: pi
    };
}


/*
 * Mean.
 */
function calculateMean(values) {

    if (values.length === 0) {
        return 0;
    }


    return (
        values.reduce(
            (sum, value) => sum + value,
            0
        )
        /
        values.length
    );
}


/*
 * Median.
 */
function calculateMedian(values) {

    if (values.length === 0) {
        return 0;
    }


    const sorted =
        [...values].sort(
            (a, b) => a - b
        );


    const middle =
        Math.floor(
            sorted.length / 2
        );


    if (sorted.length % 2) {

        return sorted[middle];
    }


    return (
        sorted[middle - 1] +
        sorted[middle]
    ) / 2;
}


/*
 * Equivalent to write_summary().
 */
function makeSweepsSummary(
    msaName,
    scores,
    windowSize,
    headers,
    sequences,
    msaLength
) {
    
    const numberOfSequences =
        sequences.length;


    const maxScore =
        Math.max(...scores);


    /*
     * Python uses 1-based location.
     */
    const indexOfMax =
        scores.indexOf(maxScore) + 1;


    const meanScore =
        calculateMean(scores);


    const medianScore =
        calculateMedian(scores);


    const maxMeanDivision =
        meanScore === 0
            ? -1
            : maxScore / meanScore;


    const maxMedianDivision =
        medianScore === 0
            ? -1
            : maxScore / medianScore;


    const relativeLocationOfPeak =
        indexOfMax / msaLength;


    const centrality =
        Math.min(
            1 - relativeLocationOfPeak,
            relativeLocationOfPeak
        ) * 2;


    const diversity =
    getAveragePairwiseDistanceAndPi(
        headers,
        sequences,
        msaLength
    );


    let above05 = 0;
    let above25 = 0;
    let above50 = 0;
    let above75 = 0;
    let above95 = 0;


    if (maxScore > 0.95) {

        above05 =
            above25 =
            above50 =
            above75 =
            above95 = 1;
    }
    else if (maxScore > 0.75) {

        above05 =
            above25 =
            above50 =
            above75 = 1;
    }
    else if (maxScore > 0.5) {

        above05 =
            above25 =
            above50 = 1;
    }
    else if (maxScore > 0.25) {

        above05 =
            above25 = 1;
    }
    else if (maxScore > 0.05) {

        above05 = 1;
    }


    const metadata = [
        msaName,
        maxScore,
        numberOfSequences,
        centrality,
        msaLength,
        windowSize,
        indexOfMax,
        meanScore,
        medianScore,
        maxMeanDivision,
        maxMedianDivision,
        relativeLocationOfPeak,
        diversity.apd,
        diversity.pi,
        above95,
        above75,
        above50,
        above25,
        above05
    ];

    const data =
    metadata[0]
    +
    ","
    +
    metadata
        .slice(1)
        .map(
            value =>
                Math.abs(value).toFixed(4)
        )
        .join(",");

    return data + "\n";
}


/*
 * Main browser equivalent of compute_sweeps_score().
 */
function computeSweepsScore(
    msaText,
    homoplasyText,
    windowSize = 50,
    msaName = "msa_fixed.fasta"
) {

    console.log(
        "Starting computeSweepsScore..."
    );


    const msa =
        parseSweepsMsa(
            msaText
        );


    const homoplasy =
        parseHomoplasy(
            homoplasyText
        );


    /*
     * Same consistency requirement as Python.
     */
    if (
        msa.msaLength !==
        homoplasy.length
    ) {

        throw new Error(
            "MSA length is inconsistent with " +
            "homoplasy length " +
            "(they should be both of the same length). " +
            "MSA: " +
            msa.msaLength +
            ", homoplasy: " +
            homoplasy.length
        );
    }


    if (
        !Number.isInteger(windowSize) ||
        windowSize < 2
    ) {

        throw new Error(
            "Window size must be an integer >= 2."
        );
    }


    const msaColumns =
        getAlignmentColumns(
            msa.sequences,
            msa.msaLength
        );


    const scores =
        getContigScores(
            msaColumns,
            msa.msaLength,
            homoplasy,
            windowSize
        );


    /*
     * Equivalent to sweeps_scores.txt.
     */
    const scoresText =
        scores
            .map(
                score =>
                    Math.abs(score)
                        .toFixed(4)
            )
            .join("\n")
        +
        "\n";


    /*
     * Equivalent to sweeps_summary.txt.
     */
    const summaryData = makeSweepsSummary(
    msaName,
    scores,
    windowSize,
    msa.headers,
    msa.sequences,
    msa.msaLength
    );
    
    const summaryHeader =
        "msa_name,max_score,number_of_sequences,centrality,msa_length," +
        "window_size,index_of_max,mean_score,median_score,max_mean_division," +
        "max_median_division,relative_location_of_peak,apd,pi," +
        "above95,above75,above50,above25,above05";
    
    const summaryText =
        summaryHeader + "\n" + summaryData;


    console.log(
        "computeSweepsScore finished."
    );

    const plotDataUrl =
        createSweepsPlot(
            scores,
            msaName,
            msa.sequences.length
        );
    
    return {
        scores,
        scoresText,
        summaryText,
        plotDataUrl
    };
}

function createSweepsPlot(scores, msaName, numberOfSequences) {

    const canvas = document.createElement("canvas");
    canvas.width = 1000;
    canvas.height = 600;

    const ctx = canvas.getContext("2d");

    const left = 80;
    const right = 30;
    const top = 70;
    const bottom = 70;

    const width = canvas.width - left - right;
    const height = canvas.height - top - bottom;

    const maxScore = 1; // Math.max(...scores);
    const minScore = 0; // Math.min(...scores);

    // White background
    ctx.fillStyle = "white";
    ctx.fillRect(0, 0, canvas.width, canvas.height);

    // Title
    ctx.fillStyle = "black";
    ctx.font = "20px Arial";
    ctx.textAlign = "center";

    ctx.fillText(
        `SINCOPA Scores for ${msaName} (across ${numberOfSequences} sequences)`,
        canvas.width / 2,
        30
    );

    // Axes
    ctx.strokeStyle = "black";
    ctx.lineWidth = 1;

    ctx.beginPath();
    ctx.moveTo(left, top);
    ctx.lineTo(left, top + height);
    ctx.lineTo(left + width, top + height);
    ctx.stroke();

    // Score line
    ctx.strokeStyle = "blue";
    ctx.lineWidth = 1;

    ctx.beginPath();

    scores.forEach((score, i) => {

        const x =
            left +
            (i / (scores.length - 1)) * width;

        const y =
            top +
            height -
            ((score - minScore) /
                (maxScore - minScore || 1)) *
                height;

        if (i === 0) {
            ctx.moveTo(x, y);
        } else {
            ctx.lineTo(x, y);
        }
    });

    ctx.stroke();

    // X axis label
    ctx.fillStyle = "black";
    ctx.font = "16px Arial";
    ctx.textAlign = "center";

    ctx.fillText(
        "Window #",
        left + width / 2,
        canvas.height - 20
    );

    // Y axis label
    ctx.save();

    ctx.translate(
        22,
        top + height / 2
    );

    ctx.rotate(-Math.PI / 2);

    ctx.textAlign = "center";

    ctx.fillText(
        "SINCOPA Score",
        0,
        0
    );

    ctx.restore();

    // X axis tick values
    ctx.font = "12px Arial";

    for (let i = 0; i <= 5; i++) {

        const index =
            Math.round(
                i * (scores.length - 1) / 5
            );

        const x =
            left + i * width / 5;

        ctx.textAlign = "center";

        ctx.fillText(
            index.toString(),
            x,
            top + height + 20
        );
    }

    // Y axis tick values
    for (let i = 0; i <= 5; i++) {

        const value =
            minScore +
            i * (maxScore - minScore) / 5;

        const y =
            top +
            height -
            i * height / 5;

        ctx.textAlign = "right";

        ctx.fillText(
            value.toFixed(2),
            left - 8,
            y + 4
        );
    }

    return canvas.toDataURL("image/png");
}
