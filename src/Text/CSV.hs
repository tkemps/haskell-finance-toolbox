{-# LANGUAGE FlexibleInstances,TypeSynonymInstances #-}

-----------------------------------------------------------------------------
-- |
-- Module      :  CRP.CSV
-- Copyright   :  (c) SKS Unternehmensberatung GmbH & Co. KG, 2012
-- License     :  see the file libraries/base/LICENSE
-- 
-- Maintainer  :  christian.haggert@sks-ub.de
-- stability   :  experimental
-- Portability :  portable
--
-- The module CSV implements functions based on parsec that allow to parse
-- a csv (comma seperated values) file. The seperators have to be commas ","
-- and decimal numbers must use decimal points (no commas like in e.g. German).
--
-----------------------------------------------------------------------------

module Text.CSV (parseCSV,csv,CSV) where

import Data.List (intercalate)
import Text.Parsec

csvFile = endBy line eol
line = sepBy cell (char ',')
cell = quotedCell <|> many (noneOf ",\n\r")

quotedCell = do 
  char '"'
  content <- many quotedChar
  char '"' <?> "quote at end of cell"
  return content

quotedChar = noneOf "\"" <|> try (string "\"\"" >> return '"')

eol = try (string "\n\r")
      <|> try (string "\r\n")
      <|> string "\n"
      <|> string "\r"
      <?> "end of line"

parseCSV :: String -> String -> Either ParseError [[String]]
parseCSV file input = parse csvFile file input

-- *Umwandlung in CSV

-- |Wandle eine Liste von Überschriften und eine Liste von Tupeln in das CSV-Format 
-- um. Falls die Liste der Überschriften leer ist, wird die ansonsten erste Zeile 
-- für die Überschriften weggelassen.
csv :: CSVTuple a => [String] -> [a] -> String
csv headers lines = intercalate "\n" allLines
    where h = intercalate "\t" headers
          b = map (\tp -> intercalate "\t" (tplToStrs tp)) lines
          allLines  = if null headers then b else h:b

class CSV a where
    csvToStr :: a -> String
instance CSV Char where
    csvToStr c = '"':c:"\""
instance CSV String where
    csvToStr s = '"':s++"\""
instance CSV Int where
    csvToStr = show
instance CSV Double where
    csvToStr = show

class CSVTuple a where
    tplToStrs :: a -> [String]
instance (CSV a,CSV b) => CSVTuple (a,b) where
    tplToStrs (a,b) = [csvToStr a,csvToStr b]
instance (CSV a,CSV b,CSV c) => CSVTuple (a,b,c) where
    tplToStrs (a,b,c) = [csvToStr a,csvToStr b,csvToStr c]
instance (CSV a,CSV b,CSV c,CSV d) => CSVTuple (a,b,c,d) where
    tplToStrs (a,b,c,d) = [csvToStr a,csvToStr b,csvToStr c,csvToStr d]
instance (CSV a,CSV b,CSV c,CSV d,CSV e) => CSVTuple (a,b,c,d,e) where
    tplToStrs (a,b,c,d,e) = [csvToStr a,csvToStr b,csvToStr c,csvToStr d,csvToStr e]
instance (CSV a,CSV b,CSV c,CSV d,CSV e,CSV f) => CSVTuple (a,b,c,d,e,f) where
    tplToStrs (a,b,c,d,e,f) = [csvToStr a,csvToStr b,csvToStr c,csvToStr d,csvToStr e,csvToStr f]

-- *Tools

-- |Ersetze jedes Vorkommen von @x@ in der Liste @as@ durch @y@.
replace x y as = map (\a -> if a==x then y else a) as

-- |Umwandlung einer ganzen Zahl in eine Fließkommazahl
dbl :: Int -> Double
dbl = fromRational . toRational

-- |Wandle Gelitkommazahlen in einen String im deutschen Format
showG :: Double -> String
showG x = replace '.' ',' (show x)

tp1lst (a) = [a]
tp2lst (a,b) = [a,b]
tp3lst (a,b,c) = [a,b,c]
tp4lst (a,b,c,d) = [a,b,c,d]
tp5lst (a,b,c,d,e) = [a,b,c,d,e]
tp6lst (a,b,c,d,e,f) = [a,b,c,d,e,f]
tp7lst (a,b,c,d,e,f,g) = [a,b,c,d,e,f,g]
