/*
 * Tong Shu Li
 * 2015-07-02
 *
 * This program uses LingPipe 4.1.0 to split an abstract
 * down to individual sentences.
 *
 * It reads input from a file and outputs the sentences
 * to a file, one sentence per line.
 *
 * Input and output files are passed to the program as
 * command line arguments.
 */
import com.aliasi.chunk.Chunk;
import com.aliasi.chunk.Chunking;

import com.aliasi.sentences.MedlineSentenceModel;
import com.aliasi.sentences.SentenceChunker;
import com.aliasi.sentences.SentenceModel;

import com.aliasi.tokenizer.IndoEuropeanTokenizerFactory;
import com.aliasi.tokenizer.TokenizerFactory;

import com.aliasi.util.Files;

// --------------- end LingPipe libraries -------------

import java.io.*;

import java.util.Iterator;
import java.util.Set;

public class SplitAbstract
{
    static final TokenizerFactory TOKENIZER_FACTORY = IndoEuropeanTokenizerFactory.INSTANCE;
    static final SentenceModel SENTENCE_MODEL = new MedlineSentenceModel();
    static final SentenceChunker SENTENCE_CHUNKER = new SentenceChunker(TOKENIZER_FACTORY, SENTENCE_MODEL);

    public static void main(String[] args) throws IOException
    {
        File in_file = new File(args[0]);
        String text = Files.readFromFile(in_file, "UTF-8");

        Chunking chunking = SENTENCE_CHUNKER.chunk(text.toCharArray(), 0, text.length());

        Set<Chunk> sentences = chunking.chunkSet();

        try
        {
            File out_file = new File(args[1]);
            if (!out_file.exists())
                out_file.createNewFile();

            FileWriter writer = new FileWriter(out_file);

            Iterator<Chunk> it = sentences.iterator();
            while (it.hasNext())
            {
                Chunk sentence = it.next();
                int start = sentence.start();
                int end = sentence.end();

                writer.write(text.substring(start, end));

                writer.write("\n");
                writer.flush();
            }

            writer.close();
        }
        catch (IOException ex)
        {
            System.out.println("Could not write to file");
        }
    }
}
