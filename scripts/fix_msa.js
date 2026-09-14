function parseFasta(msaText) {
    const sequences = {};
    const order = [];
    let currentName = null;

    for (const rawLine of msaText.split(/\r?\n/)) {
        const line = rawLine.trim();

        if (!line) continue;

        if (line.startsWith(">")) {
            currentName = line.substring(1).trim().split(/\s+/)[0];

            if (!currentName)
                throw new Error("Illegal FASTA header.");

            if (Object.prototype.hasOwnProperty.call(sequences, currentName))
                throw new Error("Duplicate FASTA sequence name: " + currentName);

            sequences[currentName] = "";
            order.push(currentName);
        }
        else {
            if (currentName === null)
                throw new Error("Illegal FASTA: sequence found before first header.");

            sequences[currentName] += line.toUpperCase();
        }
    }

    if (order.length === 0)
        throw new Error("No sequences found in MSA.");

    return { sequences, order };
}


function getFirstColumnWithoutGap(sequences, order, indexes) {
    for (const col of indexes) {
        let gapFound = false;

        for (const strain of order) {
            if (sequences[strain][col] === "-") {
                gapFound = true;
                break;
            }
        }

        if (!gapFound)
            return col;
    }

    return null;
}


function fixMsa(msaText, minimalLength = 300) {
    const { sequences, order } = parseFasta(msaText);

    const msaLength = sequences[order[0]].length;

    for (const strain of order) {
        if (sequences[strain].length !== msaLength) {
            throw new Error(
                "Illegal MSA. Not all sequences are of the same length. " +
                strain + " sequence length is " +
                sequences[strain].length +
                " where others are of length " +
                msaLength + "."
            );
        }
    }

    const leftIndexes =
        Array.from({ length: msaLength }, (_, i) => i);

    const rightIndexes =
        Array.from({ length: msaLength }, (_, i) => msaLength - 1 - i);

    const firstLeft =
        getFirstColumnWithoutGap(sequences, order, leftIndexes);

    const firstRight =
        getFirstColumnWithoutGap(sequences, order, rightIndexes);

    if (firstLeft === null || firstRight === null || firstLeft > firstRight)
        throw new Error("Could not find MSA boundaries without gaps.");

    for (const strain of order) {
        sequences[strain] =
            sequences[strain].slice(firstLeft, firstRight + 1);
    }

    const trimmedLength = firstRight - firstLeft + 1;

    if (trimmedLength < minimalLength) {
        throw new Error(
            "MSA is too short after trimming. It is " +
            trimmedLength + " bp; minimum is " +
            minimalLength + " bp."
        );
    }

    const legalChars = new Set(["A", "C", "G", "T", "-"]);

    for (let col = 0; col < trimmedLength; col++) {
        const counts = {};

        for (const strain of order) {
            const ch = sequences[strain][col];

            if (legalChars.has(ch))
                counts[ch] = (counts[ch] || 0) + 1;
        }

        if (Object.keys(counts).length === 0)
            throw new Error("Column " + (col + 1) + " contains no legal characters.");

        let majorAllele = null;
        let maxCount = -1;

        for (const ch of ["A", "C", "G", "T", "-"]) {
            const count = counts[ch] || 0;

            if (count > maxCount) {
                maxCount = count;
                majorAllele = ch;
            }
        }

        for (const strain of order) {
            if (!legalChars.has(sequences[strain][col])) {
                sequences[strain] =
                    sequences[strain].substring(0, col) +
                    majorAllele +
                    sequences[strain].substring(col + 1);
            }
        }
    }

    let fixedAlignment = "";

    for (const strain of order) {
        fixedAlignment +=
            ">" + strain + "\n" +
            sequences[strain] + "\n";
    }

    return {
        fixedMsa: fixedAlignment,
        originalLength: msaLength,
        trimmedLength: trimmedLength
    };
}
