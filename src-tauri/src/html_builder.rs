use std::fmt::Display;

#[derive(Debug, Hash, Clone)]
pub struct HtmlElement {
    tag: HtmlTag,
    header: Vec<(String, Option<String>)>,
    content: Vec<HtmlContent>,
}

impl HtmlElement {
    pub fn new(tag: HtmlTag) -> Self {
        HtmlElement {
            tag,
            header: Vec::new(),
            content: Vec::new(),
        }
    }

    pub fn class(mut self, classes: impl Into<String>) -> Self {
        self.header
            .push(("class".to_string(), Some(classes.into())));
        self
    }

    pub fn id(mut self, id: impl Into<String>) -> Self {
        self.header.push(("id".to_string(), Some(id.into())));
        self
    }

    pub fn header(mut self, title: impl Into<String>, value: impl Into<String>) -> Self {
        self.header
            .push((title.into(), Some(value.into().replace('\'', "\""))));
        self
    }

    pub fn header2(mut self, title: impl Into<String>) -> Self {
        self.header.push((title.into(), None));
        self
    }

    pub fn content(mut self, content: impl Into<HtmlContent>) -> Self {
        content.into().add_to(&mut self.content);
        self
    }
}

impl<A: Into<HtmlContent>> Extend<A> for HtmlElement {
    fn extend<T: IntoIterator<Item = A>>(&mut self, iter: T) {
        self.content.extend(iter.into_iter().map(|v| v.into()))
    }
}

impl Display for HtmlElement {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "<{}{}>",
            self.tag,
            self.header.iter().fold(String::new(), |acc, header| acc
                + " "
                + &header.0
                + &header
                    .1
                    .as_ref()
                    .map_or(String::new(), |v| format!("='{v}'")))
        )
        .unwrap();
        if !self.content.is_empty() {
            for item in &self.content {
                write!(f, "{item}").unwrap();
            }
        }
        write!(f, "</{}>", self.tag).unwrap();
        Ok(())
    }
}

pub trait ToHtmlContent {
    fn add_to(self, content: &mut Vec<HtmlContent>);
}

impl ToHtmlContent for HtmlContent {
    fn add_to(self, content: &mut Vec<HtmlContent>) {
        content.push(self)
    }
}

impl<T: IntoIterator<Item = T2>, T2: Into<HtmlContent>> ToHtmlContent for T {
    fn add_to(self, content: &mut Vec<HtmlContent>) {
        content.extend(self.into_iter().map(|i| i.into()))
    }
}

//impl<T: Into<HtmlContent>> ToHtmlContent for T {
//    fn add_to(self, content: &mut Vec<HtmlContent>) {
//        content.push(self.into())
//    }
//}

#[derive(Debug, Hash, Clone)]
pub enum HtmlContent {
    Text(String),
    Html(HtmlElement),
}

impl From<String> for HtmlContent {
    fn from(value: String) -> Self {
        HtmlContent::Text(value) // TODO: escape all html thingies
    }
}

impl From<HtmlElement> for HtmlContent {
    fn from(value: HtmlElement) -> Self {
        HtmlContent::Html(value)
    }
}

impl Display for HtmlContent {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Text(t) => write!(f, "{t}"),
            Self::Html(t) => write!(f, "{t}"),
        }
    }
}

#[derive(Debug, PartialEq, Eq, Hash, Clone, Copy)]
#[allow(non_camel_case_types)]
pub enum HtmlTag {
    /// Inline text semantics - The <a> HTML element (or anchor element), with its href attribute, creates a hyperlink to web pages, files, email addresses, locations in the same page, or anything else a URL can address.
    a,
    /// Inline text semantics - The <abbr> HTML element represents an abbreviation or acronym.
    abbr,
    /// Obsolete and deprecated - The <acronym> HTML element allows authors to clearly indicate a sequence of characters that compose an acronym or abbreviation for a word.
    acronym,
    /// Content sectioning - The <address> HTML element indicates that the enclosed HTML provides contact information for a person or people, or for an organization.
    address,
    /// Obsolete and deprecated - The obsolete HTML Applet Element (<applet>) embeds a Java applet into the document; this element has been deprecated in favor of object.
    applet,
    /// Image and multimedia - The <area> HTML element defines an area inside an image map that has predefined clickable areas. An image map allows geometric areas on an image to be associated with Hyperlink.
    area,
    /// Content sectioning - The <article> HTML element represents a self-contained composition in a document, page, application, or site, which is intended to be independently distributable or reusable (e.g., in syndication). Examples include: a forum post, a magazine or newspaper article, or a blog entry, a product card, a user-submitted comment, an interactive widget or gadget, or any other independent item of content.
    article,
    /// Content sectioning - The <aside> HTML element represents a portion of a document whose content is only indirectly related to the document's main content. Asides are frequently presented as sidebars or call-out boxes.
    aside,
    /// Image and multimedia - The <audio> HTML element is used to embed sound content in documents. It may contain one or more audio sources, represented using the src attribute or the source element: the browser will choose the most suitable one. It can also be the destination for streamed media, using a MediaStream.
    audio,
    /// Inline text semantics - The <b> HTML element is used to draw the reader's attention to the element's contents, which are not otherwise granted special importance. This was formerly known as the Boldface element, and most browsers still draw the text in boldface. However, you should not use <b> for styling text; instead, you should use the CSS font-weight property to create boldface text, or the strong element to indicate that text is of special importance.
    b,
    /// Inline text semantics - The <bdi> HTML element tells the browser's bidirectional algorithm to treat the text it contains in isolation from its surrounding text. It's particularly useful when a website dynamically inserts some text and doesn't know the directionality of the text being inserted.
    bdi,
    /// Inline text semantics - The <bdo> HTML element overrides the current directionality of text, so that the text within is rendered in a different direction.
    bdo,
    /// Obsolete and deprecated - The <bgsound> HTML element is deprecated. It sets up a sound file to play in the background while the page is used; use audio instead.
    bgsound,
    /// Obsolete and deprecated - The <big> HTML deprecated element renders the enclosed text at a font size one level larger than the surrounding text (medium becomes large, for example). The size is capped at the browser's maximum permitted font size.
    big,
    /// Obsolete and deprecated - The <blink> HTML element is a non-standard element which causes the enclosed text to flash slowly.
    blink,
    /// Text content - The <blockquote> HTML element indicates that the enclosed text is an extended quotation. Usually, this is rendered visually by indentation (see Notes for how to change it). A URL for the source of the quotation may be given using the cite attribute, while a text representation of the source can be given using the cite element.
    blockquote,
    /// Sectioning root - The <body> HTML element represents the content of an HTML document. There can be only one <body> element in a document.
    body,
    /// Inline text semantics - The <br> HTML element produces a line break in text (carriage-return). It is useful for writing a poem or an address, where the division of lines is significant.
    br,
    /// Forms - The <button> HTML element is an interactive element activated by a user with a mouse, keyboard, finger, voice command, or other assistive technology. Once activated, it then performs a programmable action, such as submitting a form or opening a dialog.
    button,
    /// Scripting - Use the HTML <canvas> element with either the canvas scripting API or the WebGL API to draw graphics and animations.
    canvas,
    /// Table content - The <caption> HTML element specifies the caption (or title) of a table.
    caption,
    /// Obsolete and deprecated - The <center> HTML element is a block-level element that displays its block-level or inline contents centered horizontally within its containing element. The container is usually, but isn't required to be, body.
    center,
    /// Inline text semantics - The <cite> HTML element is used to describe a reference to a cited creative work, and must include the title of that work. The reference may be in an abbreviated form according to context-appropriate conventions related to citation metadata.
    cite,
    /// Inline text semantics - The <code> HTML element displays its contents styled in a fashion intended to indicate that the text is a short fragment of computer code. By default, the content text is displayed using the user agent default monospace font.
    code,
    /// Table content - The <col> HTML element defines a column within a table and is used for defining common semantics on all common cells. It is generally found within a colgroup element.
    col,
    /// Table content - The <colgroup> HTML element defines a group of columns within a table.
    colgroup,
    /// Obsolete and deprecated - The <content> HTML element—an obsolete part of the Web Components suite of technologies—was used inside of Shadow DOM as an insertion point, and wasn't meant to be used in ordinary HTML. It has now been replaced by the slot element, which creates a point in the DOM at which a shadow DOM can be inserted.
    content,
    /// Inline text semantics - The <data> HTML element links a given piece of content with a machine-readable translation. If the content is time- or date-related, the time element must be used.
    data,
    /// Forms - The <datalist> HTML element contains a set of option elements that represent the permissible or recommended options available to choose from within other controls.
    datalist,
    /// Text content - The <dd> HTML element provides the description, definition, or value for the preceding term (dt) in a description list (dl).
    dd,
    /// Demarcating edits - The <del> HTML element represents a range of text that has been deleted from a document. This can be used when rendering "track changes" or source code diff information, for example. The ins element can be used for the opposite purpose: to indicate text that has been added to the document.
    del,
    /// Interactive elements - The <details> HTML element creates a disclosure widget in which information is visible only when the widget is toggled into an "open" state. A summary or label must be provided using the summary element.
    details,
    /// Inline text semantics - The <dfn> HTML element is used to indicate the term being defined within the context of a definition phrase or sentence. The p element, the dt/dd pairing, or the section element which is the nearest ancestor of the <dfn> is considered to be the definition of the term.
    dfn,
    /// Interactive elements - The <dialog> HTML element represents a dialog box or other interactive component, such as a dismissible alert, inspector, or subwindow.
    dialog,
    /// Obsolete and deprecated - The <dir> HTML element is used as a container for a directory of files and/or folders, potentially with styles and icons applied by the user agent. Do not use this obsolete element; instead, you should use the ul element for lists, including lists of files.
    dir,
    /// Text content - The <div> HTML element is the generic container for flow content. It has no effect on the content or layout until styled in some way using CSS (e.g. styling is directly applied to it, or some kind of layout model like Flexbox is applied to its parent element).
    div,
    /// Text content - The <dl> HTML element represents a description list. The element encloses a list of groups of terms (specified using the dt element) and descriptions (provided by dd elements). Common uses for this element are to implement a glossary or to display metadata (a list of key-value pairs).
    dl,
    /// Text content - The <dt> HTML element specifies a term in a description or definition list, and as such must be used inside a dl element. It is usually followed by a dd element; however, multiple <dt> elements in a row indicate several terms that are all defined by the immediate next dd element.
    dt,
    /// Inline text semantics - The <em> HTML element marks text that has stress emphasis. The <em> element can be nested, with each level of nesting indicating a greater degree of emphasis.
    em,
    /// Embedded content - The <embed> HTML element embeds external content at the specified point in the document. This content is provided by an external application or other source of interactive content such as a browser plug-in.
    embed,
    /// Forms - The <fieldset> HTML element is used to group several controls as well as labels (label) within a web form.
    fieldset,
    /// Text content - The <figcaption> HTML element represents a caption or legend describing the rest of the contents of its parent figure element.
    figcaption,
    /// Text content - The <figure> HTML element represents self-contained content, potentially with an optional caption, which is specified using the figcaption element. The figure, its caption, and its contents are referenced as a single unit.
    figure,
    /// Obsolete and deprecated - The <font> HTML element defines the font size, color and face for its content.
    font,
    /// Content sectioning - The <footer> HTML element represents a footer for its nearest ancestor sectioning content or sectioning root element. A <footer> typically contains information about the author of the section, copyright data or links to related documents.
    footer,
    /// Forms - The <form> HTML element represents a document section containing interactive controls for submitting information.
    form,
    /// Obsolete and deprecated - The <frame> HTML element defines a particular area in which another HTML document can be displayed. A frame should be used within a frameset.
    frame,
    /// Obsolete and deprecated - The <frameset> HTML element is used to contain frame elements.
    frameset,
    /// Content sectioning - The <h1> to <h6> HTML elements represent six levels of section headings. <h1> is the highest section level and <h6> is the lowest.
    h1,
    /// Content sectioning - The <h1> to <h6> HTML elements represent six levels of section headings. <h1> is the highest section level and <h6> is the lowest.
    h2,
    /// Content sectioning - The <h1> to <h6> HTML elements represent six levels of section headings. <h1> is the highest section level and <h6> is the lowest.
    h3,
    /// Content sectioning - The <h1> to <h6> HTML elements represent six levels of section headings. <h1> is the highest section level and <h6> is the lowest.
    h4,
    /// Content sectioning - The <h1> to <h6> HTML elements represent six levels of section headings. <h1> is the highest section level and <h6> is the lowest.
    h5,
    /// Content sectioning - The <h1> to <h6> HTML elements represent six levels of section headings. <h1> is the highest section level and <h6> is the lowest.
    h6,
    /// Document metadata - The <head> HTML element contains machine-readable information (metadata) about the document, like its title, scripts, and style sheets.
    head,
    /// Content sectioning - The <header> HTML element represents introductory content, typically a group of introductory or navigational aids. It may contain some heading elements but also a logo, a search form, an author name, and other elements.
    header,
    /// Text content - The <hr> HTML element represents a thematic break between paragraph-level elements: for example, a change of scene in a story, or a shift of topic within a section.
    hr,
    /// Main - The <html> HTML element represents the root (top-level element) of an HTML document, so it is also referred to as the root element. All other elements must be descendants of this element.
    html,
    /// Inline text semantics - The <i> HTML element represents a range of text that is set off from the normal text for some reason, such as idiomatic text, technical terms, taxonomical designations, among others. Historically, these have been presented using italicized type, which is the original source of the <i> naming of this element.
    i,
    /// Embedded content - The <iframe> HTML element represents a nested browsing context, embedding another HTML page into the current one.
    iframe,
    /// Obsolete and deprecated - The <image> HTML element is an ancient and poorly supported precursor to the img element. It should not be used.
    image,
    /// Image and multimedia - The <img> HTML element embeds an image into the document.
    img,
    /// Forms - The <input> HTML element is used to create interactive controls for web-based forms in order to accept data from the user; a wide variety of types of input data and control widgets are available, depending on the device and user agent. The <input> element is one of the most powerful and complex in all of HTML due to the sheer number of combinations of input types and attributes.
    input,
    /// Demarcating edits - The <ins> HTML element represents a range of text that has been added to a document. You can use the del element to similarly represent a range of text that has been deleted from the document.
    ins,
    /// Inline text semantics - The <kbd> HTML element represents a span of inline text denoting textual user input from a keyboard, voice input, or any other text entry device. By convention, the user agent defaults to rendering the contents of a <kbd> element using its default monospace font, although this is not mandated by the HTML standard.
    kbd,
    /// Obsolete and deprecated - The <keygen> HTML element exists to facilitate generation of key material, and submission of the public key as part of an HTML form. This mechanism is designed for use with Web-based certificate management systems. It is expected that the <keygen> element will be used in an HTML form along with other information needed to construct a certificate request, and that the result of the process will be a signed certificate.
    keygen,
    /// Forms - The <label> HTML element represents a caption for an item in a user interface.
    label,
    /// Forms - The <legend> HTML element represents a caption for the content of its parent fieldset.
    legend,
    /// Text content - The <li> HTML element is used to represent an item in a list. It must be contained in a parent element: an ordered list (ol), an unordered list (ul), or a menu (menu). In menus and unordered lists, list items are usually displayed using bullet points. In ordered lists, they are usually displayed with an ascending counter on the left, such as a number or letter.
    li,
    /// Document metadata - The <link> HTML element specifies relationships between the current document and an external resource. This element is most commonly used to link to CSS, but is also used to establish site icons (both "favicon" style icons and icons for the home screen and apps on mobile devices) among other things.
    link,
    /// Content sectioning - The <main> HTML element represents the dominant content of the body of a document. The main content area consists of content that is directly related to or expands upon the central topic of a document, or the central functionality of an application.
    main,
    /// Image and multimedia - The <map> HTML element is used with area elements to define an image map (a clickable link area).
    map,
    /// Inline text semantics - The <mark> HTML element represents text which is marked or highlighted for reference or notation purposes, due to the marked passage's relevance or importance in the enclosing context.
    mark,
    /// Obsolete and deprecated - The <marquee> HTML element is used to insert a scrolling area of text. You can control what happens when the text reaches the edges of its content area using its attributes.
    marquee,
    /// The top-level element in MathML is <math>. Every valid MathML instance must be wrapped in <math> tags. In addition you must not nest a second <math> element in another, but you can have an arbitrary number of other child elements in it.
    math,
    /// Text content - The <menu> HTML element is described in the HTML specification as a semantic alternative to ul, but treated by browsers (and exposed through the accessibility tree) as no different than ul. It represents an unordered list of items (which are represented by li elements).
    menu,
    /// Obsolete and deprecated - The <menuitem> HTML element represents a command that a user is able to invoke through a popup menu. This includes context menus, as well as menus that might be attached to a menu button.
    menuitem,
    /// Document metadata - The <meta> HTML element represents Metadata that cannot be represented by other HTML meta-related elements, like base, link, script, style or title.
    meta,
    /// Forms - The <meter> HTML element represents either a scalar value within a known range or a fractional value.
    meter,
    /// Content sectioning - The <nav> HTML element represents a section of a page whose purpose is to provide navigation links, either within the current document or to other documents. Common examples of navigation sections are menus, tables of contents, and indexes.
    nav,
    /// Obsolete and deprecated - The <nobr> HTML element prevents the text it contains from automatically wrapping across multiple lines, potentially resulting in the user having to scroll horizontally to see the entire width of the text.
    nobr,
    /// Obsolete and deprecated - The <noembed> HTML element is an obsolete, non-standard way to provide alternative, or "fallback", content for browsers that do not support the embed element or do not support the type of embedded content an author wishes to use. This element was deprecated in HTML 4.01 and above in favor of placing fallback content between the opening and closing tags of an object element.
    noembed,
    /// Obsolete and deprecated - The <noframes> HTML element provides content to be presented in browsers that don't support (or have disabled support for) the frame element. Although most commonly-used browsers support frames, there are exceptions, including certain special-use browsers including some mobile browsers, as well as text-mode browsers.
    noframes,
    /// Scripting - The <noscript> HTML element defines a section of HTML to be inserted if a script type on the page is unsupported or if scripting is currently turned off in the browser.
    noscript,
    /// Text content - The <ol> HTML element represents an ordered list of items — typically rendered as a numbered list.
    ol,
    /// Forms - The <optgroup> HTML element creates a grouping of options within a select element.
    optgroup,
    /// Forms - The <option> HTML element is used to define an item contained in a select, an optgroup, or a datalist element. As such, <option> can represent menu items in popups and other lists of items in an HTML document.
    option,
    /// Forms - The <output> HTML element is a container element into which a site or app can inject the results of a calculation or the outcome of a user action.
    output,
    /// Text content - The <p> HTML element represents a paragraph. Paragraphs are usually represented in visual media as blocks of text separated from adjacent blocks by blank lines and/or first-line indentation, but HTML paragraphs can be any structural grouping of related content, such as images or form fields.
    p,
    /// Obsolete and deprecated - The <param> HTML element defines parameters for an object element.
    param,
    /// Embedded content - The <picture> HTML element contains zero or more source elements and one img element to offer alternative versions of an image for different display/device scenarios.
    picture,
    /// Obsolete and deprecated - The <plaintext> HTML element renders everything following the start tag as raw text, ignoring any following HTML. There is no closing tag, since everything after it is considered raw text.
    plaintext,
    /// Embedded content - The <portal> HTML element enables the embedding of another HTML page into the current one for the purposes of allowing smoother navigation into new pages.
    portal,
    /// Text content - The <pre> HTML element represents preformatted text which is to be presented exactly as written in the HTML file. The text is typically rendered using a non-proportional, or monospaced, font. Whitespace inside this element is displayed as written.
    pre,
    /// Forms - The <progress> HTML element displays an indicator showing the completion progress of a task, typically displayed as a progress bar.
    progress,
    /// Inline text semantics - The <q> HTML element indicates that the enclosed text is a short inline quotation. Most modern browsers implement this by surrounding the text in quotation marks. This element is intended for short quotations that don't require paragraph breaks; for long quotations use the blockquote element.
    q,
    /// Obsolete and deprecated - The <rb> HTML element is used to delimit the base text component of a ruby annotation, i.e. the text that is being annotated. One <rb> element should wrap each separate atomic segment of the base text.
    rb,
    /// Inline text semantics - The <rp> HTML element is used to provide fall-back parentheses for browsers that do not support display of ruby annotations using the ruby element. One <rp> element should enclose each of the opening and closing parentheses that wrap the rt element that contains the annotation's text.
    rp,
    /// Inline text semantics - The <rt> HTML element specifies the ruby text component of a ruby annotation, which is used to provide pronunciation, translation, or transliteration information for East Asian typography. The <rt> element must always be contained within a ruby element.
    rt,
    /// Obsolete and deprecated - The <rtc> HTML element embraces semantic annotations of characters presented in a ruby of rb elements used inside of ruby element. rb elements can have both pronunciation (rt) and semantic (rtc) annotations.
    rtc,
    /// Inline text semantics - The <ruby> HTML element represents small annotations that are rendered above, below, or next to base text, usually used for showing the pronunciation of East Asian characters. It can also be used for annotating other kinds of text, but this usage is less common.
    ruby,
    /// Inline text semantics - The <s> HTML element renders text with a strikethrough, or a line through it. Use the <s> element to represent things that are no longer relevant or no longer accurate. However, <s> is not appropriate when indicating document edits; for that, use the del and ins elements, as appropriate.
    s,
    /// Inline text semantics - The <samp> HTML element is used to enclose inline text which represents sample (or quoted) output from a computer program. Its contents are typically rendered using the browser's default monospaced font (such as Courier or Lucida Console).
    samp,
    /// Scripting - The <script> HTML element is used to embed executable code or data; this is typically used to embed or refer to JavaScript code. The <script> element can also be used with other languages, such as WebGL's GLSL shader programming language and JSON.
    script,
    /// Content sectioning - The <section> HTML element represents a generic standalone section of a document, which doesn't have a more specific semantic element to represent it. Sections should always have a heading, with very few exceptions.
    section,
    /// Forms - The <select> HTML element represents a control that provides a menu of options.
    select,
    /// Obsolete and deprecated - The <shadow> HTML element—an obsolete part of the Web Components technology suite—was intended to be used as a shadow DOM insertion point. You might have used it if you have created multiple shadow roots under a shadow host. It is not useful in ordinary HTML.
    shadow,
    /// Web components - The <slot> HTML element—part of the Web Components technology suite—is a placeholder inside a web component that you can fill with your own markup, which lets you create separate DOM trees and present them together.
    slot,
    /// Inline text semantics - The <small> HTML element represents side-comments and small print, like copyright and legal text, independent of its styled presentation. By default, it renders text within it one font-size smaller, such as from small to x-small.
    small,
    /// Embedded content - The <source> HTML element specifies multiple media resources for the picture, the audio element, or the video element. It is an empty element, meaning that it has no content and does not have a closing tag. It is commonly used to offer the same media content in multiple file formats in order to provide compatibility with a broad range of browsers given their differing support for image file formats and media file formats.
    source,
    /// Obsolete and deprecated - The <spacer> HTML element is an obsolete HTML element which allowed insertion of empty spaces on pages. It was devised by Netscape to accomplish the same effect as a single-pixel layout image, which was something web designers used to use to add white spaces to web pages without actually using an image. However, <spacer> no longer supported by any major browser and the same effects can now be achieved using simple CSS.
    spacer,
    /// Inline text semantics - The <span> HTML element is a generic inline container for phrasing content, which does not inherently represent anything. It can be used to group elements for styling purposes (using the class or id attributes), or because they share attribute values, such as lang. It should be used only when no other semantic element is appropriate. <span> is very much like a div element, but div is a block-level element whereas a <span> is an inline element.
    span,
    /// Obsolete and deprecated - The <strike> HTML element places a strikethrough (horizontal line) over text.
    strike,
    /// Inline text semantics - The <strong> HTML element indicates that its contents have strong importance, seriousness, or urgency. Browsers typically render the contents in bold type.
    strong,
    /// Document metadata - The <style> HTML element contains style information for a document, or part of a document. It contains CSS, which is applied to the contents of the document containing the <style> element.
    style,
    /// Inline text semantics - The <sub> HTML element specifies inline text which should be displayed as subscript for solely typographical reasons. Subscripts are typically rendered with a lowered baseline using smaller text.
    sub,
    /// Interactive elements - The <summary> HTML element specifies a summary, caption, or legend for a details element's disclosure box. Clicking the <summary> element toggles the state of the parent <details> element open and closed.
    summary,
    /// Inline text semantics - The <sup> HTML element specifies inline text which is to be displayed as superscript for solely typographical reasons. Superscripts are usually rendered with a raised baseline using smaller text.
    sup,
    /// The svg element is a container that defines a new coordinate system and viewport. It is used as the outermost element of SVG documents, but it can also be used to embed an SVG fragment inside an SVG or HTML document.
    svg,
    /// Table content - The <table> HTML element represents tabular data — that is, information presented in a two-dimensional table comprised of rows and columns of cells containing data.
    table,
    /// Table content - The <tbody> HTML element encapsulates a set of table rows (tr elements), indicating that they comprise the body of the table (table).
    tbody,
    /// Table content - The <td> HTML element defines a cell of a table that contains data. It participates in the table model.
    td,
    /// Web components - The <template> HTML element is a mechanism for holding HTML that is not to be rendered immediately when a page is loaded but may be instantiated subsequently during runtime using JavaScript.
    template,
    /// Forms - The <textarea> HTML element represents a multi-line plain-text editing control, useful when you want to allow users to enter a sizeable amount of free-form text, for example a comment on a review or feedback form.
    textarea,
    /// Table content - The <tfoot> HTML element defines a set of rows summarizing the columns of the table.
    tfoot,
    /// Table content - The <th> HTML element defines a cell as header of a group of table cells. The exact nature of this group is defined by the scope and headers attributes.
    th,
    /// Table content - The <thead> HTML element defines a set of rows defining the head of the columns of the table.
    thead,
    /// Inline text semantics - The <time> HTML element represents a specific period in time. It may include the datetime attribute to translate dates into machine-readable format, allowing for better search engine results or custom features such as reminders.
    time,
    /// Document metadata - The <title> HTML element defines the document's title that is shown in a Browser's title bar or a page's tab. It only contains text; tags within the element are ignored.
    title,
    /// Table content - The <tr> HTML element defines a row of cells in a table. The row's cells can then be established using a mix of td (data cell) and th (header cell) elements.
    tr,
    /// Image and multimedia - The <track> HTML element is used as a child of the media elements, audio and video. It lets you specify timed text tracks (or time-based data), for example to automatically handle subtitles. The tracks are formatted in WebVTT format (.vtt files) — Web Video Text Tracks.
    track,
    /// Obsolete and deprecated - The <tt> HTML element creates inline text which is presented using the user agent default monospace font face. This element was created for the purpose of rendering text as it would be displayed on a fixed-width display such as a teletype, text-only screen, or line printer.
    tt,
    /// Inline text semantics - The <u> HTML element represents a span of inline text which should be rendered in a way that indicates that it has a non-textual annotation. This is rendered by default as a simple solid underline, but may be altered using CSS.
    u,
    /// Text content - The <ul> HTML element represents an unordered list of items, typically rendered as a bulleted list.
    ul,
    /// Inline text semantics - The <var> HTML element represents the name of a variable in a mathematical expression or a programming context. It's typically presented using an italicized version of the current typeface, although that behavior is browser-dependent.
    var,
    /// Image and multimedia - The <video> HTML element embeds a media player which supports video playback into the document. You can use <video> for audio content as well, but the audio element may provide a more appropriate user experience.
    video,
    /// Inline text semantics - The <wbr> HTML element represents a word break opportunity—a position within text where the browser may optionally break a line, though its line-breaking rules would not otherwise create a break at that location.
    wbr,
    /// Obsolete and deprecated - The <xmp> HTML element renders text between the start and end tags without interpreting the HTML in between and using a monospaced font. The HTML2 specification recommended that it should be rendered wide enough to allow 80 characters per line.
    xmp,
}

impl Display for HtmlTag {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{self:?}")
    }
}
