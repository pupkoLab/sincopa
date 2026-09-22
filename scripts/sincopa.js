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
