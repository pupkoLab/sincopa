console.log("validate_input.js loaded");


function verifyFastaFormat(msaText) {

    // Same legal DNA characters accepted by the original SINCOPA
    const legalChars =
        new Set("ACGTRYSWKMBDHVNU-acgtryswkmbdhvnu-");

    const lines = msaText.split(/\r?\n/);

    if (lines.length === 0 || !lines[0].startsWith(">")) {

        const firstChar =
            lines[0] && lines[0].length > 0
                ? lines[0][0]
                : "";

        return `Illegal FASTA format. First line in MSA starts with "${firstChar}" instead of ">".`;
    }


    let previousLineWasHeader = true;
    let firstEmptyLine = null;

    for (let i = 1; i < lines.length; i++) {

        const lineNumber = i + 1;
        const line = lines[i].trim();

        // Empty line
        if (!line) {

            // Allow trailing empty lines
            const remainingNonEmpty =
                lines
                    .slice(i + 1)
                    .some(x => x.trim() !== "");

            if (remainingNonEmpty && firstEmptyLine === null) {
                firstEmptyLine = lineNumber;
            }

            continue;
        }


        // Non-empty line following an empty line
        if (firstEmptyLine !== null) {

            return `Illegal FASTA format. Line ${firstEmptyLine} in MSA is empty.`;
        }


        // Header
        if (line.startsWith(">")) {

            if (previousLineWasHeader) {

                return `Illegal FASTA format. MSA contains an empty record. Both lines ${lineNumber - 1} and ${lineNumber} start with ">".`;
            }

            previousLineWasHeader = true;
            continue;
        }


        // Sequence
        previousLineWasHeader = false;

        for (const c of line) {

            if (!legalChars.has(c)) {

                return `Illegal FASTA format. Line ${lineNumber} in MSA contains illegal DNA character "${c}".`;
            }
        }
    }


    if (previousLineWasHeader) {
        return "Illegal FASTA format. The last MSA record contains no sequence.";
    }


    return null;
}


function curateFasta(msaText) {

    // The Python version replaces "|" in headers with "_"
    const lines = msaText.split(/\r?\n/);

    return lines
        .map(line => {
            if (line.startsWith(">")) {
                return line.replace(/\|/g, "_");
            }
            return line;
        })
        .join("\n");
}


function verifyMsaIsConsistentWithTree(msaText, treeText) {

    const treeTaxa =
        new Set(getTreeLabels(treeText));

    const msaTaxa = [];

    const lines =
        msaText.split(/\r?\n/);

    for (const line of lines) {

        if (line.startsWith(">")) {

            const strain =
                line.substring(1).trim();

            msaTaxa.push(strain);
        }
    }


    for (const strain of msaTaxa) {

        if (!treeTaxa.has(strain)) {

            return `${strain} species appears in the input MSA but not in the phylogenetic tree. Please make sure the phylogenetic tree you provide contains (at least) all the species in the provided MSA.`;
        }
    }


    return null;
}


function validateInput(msaText, treeText) {

    let errorMsg;


    // FASTA validation
    errorMsg =
        verifyFastaFormat(msaText);

    if (errorMsg) {
        throw new Error(errorMsg);
    }


    // Curate MSA exactly as the old program did
    const curatedMsa =
        curateFasta(msaText);


    // Tree/MSA consistency
    errorMsg =
        verifyMsaIsConsistentWithTree(
            curatedMsa,
            treeText
        );

    if (errorMsg) {
        throw new Error(errorMsg);
    }


    return curatedMsa;
}
