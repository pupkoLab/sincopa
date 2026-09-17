console.log("fix_tree.js loaded");
let removeTaxaModule = null;

async function initializeRemoveTaxa() {
    removeTaxaModule = await createRemoveTaxaModule({
        noInitialRun: true
    });

    console.log("removeTaxa WebAssembly ready");
}

initializeRemoveTaxa();


function removeBootstrapValues(treeText) {
    // Python:
    // re.sub(r'\)\d+:', '):', tree_as_str)
    return treeText.replace(/\)\d+:/g, '):');
}


function getMsaTaxa(msaText) {
    const taxa = [];

    const regex = /^>(\S+)/gm;
    let match;

    while ((match = regex.exec(msaText)) !== null) {
        taxa.push(match[1]);
    }

    return taxa;
}


/*
 * Equivalent purpose to Bio.Phylo:
 *
 * [leaf.name for leaf in tree.get_terminals()]
 *
 * For now this extracts Newick leaf names.
 */
function getTreeLabels(treeText) {

    const labels = [];

    /*
     * Leaf labels occur following "(" or ","
     * and terminate at :, comma, or ")".
     */
    const regex = /(?:\(|,)\s*([^():,;\s]+)\s*(?=[:),])/g;

    let match;

    while ((match = regex.exec(treeText)) !== null) {
        labels.push(match[1]);
    }

    return labels;
}


async function fixTree(msaText, treeText) {

    if (!removeTaxaModule) {
        throw new Error("removeTaxa WebAssembly is still loading.");
    }

    /*
     * Step 1:
     * remove bootstrap values
     */
    const treeWithoutBootstrap =
        removeBootstrapValues(treeText);


    /*
     * Step 2:
     * get taxa from tree
     */
    const treeTaxa =
        getTreeLabels(treeWithoutBootstrap);


    /*
     * Step 3:
     * get taxa from MSA
     */
    const msaTaxa =
        getMsaTaxa(msaText);

    const msaSet =
        new Set(msaTaxa);


    /*
     * Step 4:
     * determine which tree taxa should be removed
     */
    const taxaToRemove =
        treeTaxa.filter(
            taxon => !msaSet.has(taxon)
        );

    console.log("Tree taxa:", treeTaxa);
    console.log("MSA taxa:", msaTaxa);
    console.log("Taxa to remove:", taxaToRemove);


    /*
     * Step 5:
     * put files into the Emscripten filesystem
     */
    const FS = removeTaxaModule.FS;

    const inputTreePath = "/tree.no_bootstrap";
    const taxaPath = "/taxa_to_remove.txt";
    const outputTreePath = "/tree_fixed.newick";


    /*
     * Clean old files in case user runs twice.
     */
    for (const path of [
        inputTreePath,
        taxaPath,
        outputTreePath
    ]) {
        try {
            FS.unlink(path);
        }
        catch (e) {
        }
    }


    FS.writeFile(
        inputTreePath,
        treeWithoutBootstrap
    );

    FS.writeFile(
        taxaPath,
        taxaToRemove.join("\n")
    );


    /*
     * Equivalent to:
     *
     * removeTaxa input_tree taxa_file output_tree
     */
    removeTaxaModule.callMain([
        inputTreePath,
        taxaPath,
        outputTreePath
    ]);


    /*
     * Step 6:
     * retrieve fixed tree from WASM filesystem
     */
    const fixedTree =
        FS.readFile(
            outputTreePath,
            { encoding: "utf8" }
        );


    return {
        fixedTree,
        taxaToRemove,
        treeTaxa,
        msaTaxa
    };
}


document.addEventListener(
    "DOMContentLoaded",
    function () {

        const form =
            document.getElementById("form");


        form.addEventListener(
            "submit",
            async function (event) {

                event.preventDefault();

                try {

                    const jobTitle =
                        document.getElementById("job_title").value.trim();
                    
                    const inputIsProvidedAsText =
                        document.getElementById("input_is_provided_as_text").checked;
                    
                    let msaText;
                    let treeText;
                    let msaFileName;
                    let treeFileName;
                    
                    if (inputIsProvidedAsText) {
                    
                        // Read MSA and tree directly from the textareas
                        msaText =
                            document.getElementById("msa_text").value.trim();
                    
                        treeText =
                            document.getElementById("tree_text").value.trim();
                    
                        if (!msaText) {
                            throw new Error("Please provide the multiple sequence alignment text.");
                        }
                    
                        if (!treeText) {
                            throw new Error("Please provide the phylogenetic tree text.");
                        }
                    
                        // Used on results.html and in the SINCOPA summary
                        msaFileName = "raw text";
                        treeFileName = "raw text";
                    
                    } else {
                    
                        // Read uploaded files
                        const msaFile =
                            document.getElementById("msa_file").files[0];
                    
                        const treeFile =
                            document.getElementById("tree_file").files[0];
                    
                        if (!msaFile) {
                            throw new Error("Please select a multiple sequence alignment file.");
                        }
                    
                        if (!treeFile) {
                            throw new Error("Please select a phylogenetic tree file.");
                        }
                    
                        msaText =
                            await msaFile.text();
                    
                        treeText =
                            await treeFile.text();
                    
                        msaFileName = msaFile.name;
                        treeFileName = treeFile.name;
                    }

                    /*
                     * Step 0:
                     * validate msaText, treeText
                     */
                    msaText =
                        validateInput(msaText, treeText);
                    
                    /*
                     * Step 1:
                     * fix tree using the ORIGINAL MSA
                     */
                    const treeResult =
                        await fixTree(
                            msaText,
                            treeText
                        );
                    
                    
                    /*
                     * Step 2:
                     * fix MSA
                     */
                    const msaResult =
                        fixMsa(
                            msaText,
                            300
                        );

                    /*
                     * Step 3:
                     * compute homoplasy
                     */
                    const homoplasyResult =
                        await computeHomoplasy(
                            msaResult.fixedMsa,
                            treeResult.fixedTree
                        );

                    const sweepsResult =
                    computeSweepsScore(
                        msaResult.fixedMsa,
                        homoplasyResult.homoplasyText,
                        50,
                        msaFileName
                    );
                    
                    /*
                     * Save fixed tree
                     */
                    sessionStorage.setItem(
                        "fixedTree",
                        treeResult.fixedTree
                    );
                    
                    sessionStorage.setItem(
                        "taxaToRemove",
                        JSON.stringify(
                            treeResult.taxaToRemove
                        )
                    );
                    
                    sessionStorage.setItem(
                        "treeFileName",
                        treeFileName
                    );
                    
                    
                    /*
                     * Save fixed MSA
                     */
                    sessionStorage.setItem(
                        "fixedMsa",
                        msaResult.fixedMsa
                    );
                    
                    sessionStorage.setItem(
                        "msaFileName",
                        msaFileName
                    );
                    
                    sessionStorage.setItem(
                        "msaOriginalLength",
                        msaResult.originalLength
                    );
                    
                    sessionStorage.setItem(
                        "msaTrimmedLength",
                        msaResult.trimmedLength
                    );

                    sessionStorage.setItem(
                        "inputIsProvidedAsText",
                        inputIsProvidedAsText ? "true" : "false"
                    );
                    
                    /*
                     * Save homoplasy results
                     */
                    sessionStorage.setItem(
                        "homoplasyText",
                        homoplasyResult.homoplasyText
                    );
                    
                    sessionStorage.setItem(
                        "controlText",
                        homoplasyResult.controlText
                    );

                    /*
                     * Save sweep scores results
                     */
                    
                   sessionStorage.setItem(
                        "sweepsScoresText",
                        sweepsResult.scoresText
                    );
                    
                    sessionStorage.setItem(
                        "sweepsSummaryText",
                        sweepsResult.summaryText
                    );

                    sessionStorage.setItem(
                        "sweepsPlot",
                        sweepsResult.plotDataUrl
                    );

                    sessionStorage.setItem("jobTitle", jobTitle);
                    
                    /*
                     * Go to results page
                     */
                    window.location.href =
                        "results.html";

                }
                catch (error) {

                    console.error(error);

                    alert(error.message);
                }
            }
        );
    }
);
